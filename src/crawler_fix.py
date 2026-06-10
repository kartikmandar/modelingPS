import glob
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for saving plots
from scipy.stats    import f                      as ftest
from scipy.optimize import minimize               as mini
from scipy.optimize import basinhopping           as bashop
from scipy.optimize import differential_evolution as diffevol
from scipy.optimize import shgo                   as shgo
from scipy.optimize import dual_annealing         as duann
from scipy.optimize import LinearConstraint       as LC
#from scipy.optimize import direct                 as direct

def lore(r, D, nu0, nu):
    '''
    Calculates a Lorentzian

    There is a + 1.e-20 in the denominator to
    avoid division by zero.

    Takes 

    r: ampl
    D: size
    nu0: centroid
    nu: np array where to evaluate

    Returns:

    The values
    '''

    out = (r**2 * D / np.pi) / (D**2 + (nu - nu0)**2 + 1.e-20)

    return out

class obs(object):
    '''
    Stores an observation psd and can fit it

    Has Attributes:

    freq     : the frequencies of the power spectrum
    power    : the power of the power spectrum
    err_power: the errors on the estimates
    weigh    : 1 / err^2 - used for fitting

    Internally used when searching for QPOs:

    res      : the solution without QPO
    fr       : the frequency to use for the QPO fit

    Has Methods:

    fit: main method: fits the data and QPO

    qpo_try_search: internal: fits the QPO at various frequencies
                    and keeps the one with the lowest resulting chi2

    qpo_try_max: internal: fits the QPO at the max of nu f(nu)
                 and then rifits letting everything free

    fun    : the function fitted by fit
    funQPO : the function fitted while searching for the QPO
    funQPO2: the function fitted while searching for the QPO

    make_args: internal: decides whether to pass constraints
               or bounds depending on the method of choice
               (should be a static method)

    fit_single: for testing
    '''

    def __init__(self, infile):
        '''
        infile: data file
        '''

        print('\n' + infile + '\n')

        self.infile = infile

        dat = np.genfromtxt(self.infile, skip_header=1).T

        self.freq = dat[0]
        self.power = dat[2]
        self.err_power = dat[3]
        self.weigh = 1 / self.err_power**2


    def fit_single(self, x0_given=None, plot=False):
        '''
        Fits 3 Lorentzians

        Takes:

        x0_given: an initial guess
        plot    : whether to plot the result
        '''

        chis      = []
        num       = 3
        fit_use   = None
        const_use = None
        conv      = False

        A  = np.identity(3 * num)
        lb = np.zeros(3 * num)
        ub = np.full(3 * num, 100)

        for i in range(num):
            ub[3 * i] = 1

        for i in range(1, num):
            pos = 3 * (i + 1) - 1
            A[pos, pos - 3] = -1

        const = LC(A, lb=lb, ub=ub)

        x0 = np.empty(3 * num)

        if x0_given is None:

            for i in range(num):

                ini = 3 * i

                x0[ini    ] = 0.3
                x0[ini + 1] = 0.5
                x0[ini + 2] = i * 1.e-3

        else:

            x0 = x0_given

        print(' initial conditions ', x0 )

        # min_fit = diffevol(self.fun, constraints=const,
        #                    bounds=[[0, 100] for i in range(3 * num)])
        fit_use = mini(self.fun, x0=x0, constraints=const_use,
                       method='COBYLA', options={'maxiter': 10000})

        chis.append(fit_use.fun / (len(self.freq) - 3 * num))

        print(fit_use.success)
        print('Chi2', chis[-1])

        print('results ', fit_use.x)
        

        if plot:

            fs = []
            for i in range(num):
                ini = 3 * i
                fs.append(lore(*fit_use.x[ini:ini + 3], self.freq))
                # fs.append(lore(*x0[ini:ini + 3], self.freq))
                
            fig, axi = plt.subplots(1, 1)

            axi.errorbar(self.freq,
                         self.freq * self.power,
                         yerr=self.freq * self.err_power, marker='o')

            ft = np.copy(fs[0])

            for f in fs[1:]:
                ft += f

            axi.plot(self.freq, self.freq * ft, '-')

            for f in fs:
                axi.plot(self.freq, self.freq * f, ':')

            axi.set_xscale('log')
            axi.set_yscale('log')
            axi.set_ylabel(r'$\nu f(\nu)$')
            axi.set_xlabel(r'$\nu$')
    
            fig.savefig(self.infile + '.pdf', dpi=300)
    
    def fun(self, par, full=False):
        '''
        Models n Lorentzians and subtracts them from the data

        The number of Lorentzians is deduced from the parameters.

        Takes:

        par : the parameters for the Lorentzians
        full: whether to return the residuals or the total chi2

        Returns:

        chi2 or residuals depending on full
        '''

        n = len(par) // 3

        chi = np.zeros(len(self.power))

        for i in range(n):

            ini = 3 *  i
            fin = 3 * (i + 1)

            r, D, nu0 = par[ini:fin]

            chi += lore(r, D, nu0, self.freq)

        chi = self.power - chi

        if full:
            return chi
        else:
            return np.sum(chi**2 * self.weigh)

    def funQPO(self, par, full=False):
        '''
        Models n + 1 Lorentzians and subtracts them from the data

        The number of Lorentzians is deduced from the attribute
        self.res. The frequency of the last, QPO, Lorentzian is
        taken from the attribute self.fr.

        Takes:

        par : the parameters for the Lorentzian of the QPO, ONLY
        full: whether to return the residuals or the total chi2

        Returns:

        chi2 or residuals depending on full
        '''

        n = len(self.res) // 3

        chi = np.zeros(len(self.power))

        for i in range(n):

            ini = 3 *  i
            fin = 3 * (i + 1)

            r, D, nu0 = self.res[ini:fin]

            chi += lore(r, D, nu0, self.freq)

        chi += lore(par[0], par[1], self.fr, self.freq)

        chi = self.power - chi

        if full:
            return chi
        else:
            return np.sum(chi**2 * self.weigh)

    def funQPO2(self, par, full=False):
        '''
        Models n + 1 Lorentzians and subtracts them from the data

        The number of Lorentzians is deduced from the parameters.
        The frequency of the last, QPO, Lorentzian is taken from
        the attribute self.fr.

        Takes:

        par : the parameters for the Lorentzians + the QPO
        full: whether to return the residuals or the total chi2

        Returns:

        chi2 or residuals depending on full
        '''

        n = (len(par) - 2) // 3

        chi = np.zeros(len(self.power))

        for i in range(n):

            ini = 3 *  i
            fin = 3 * (i + 1)

            r, D, nu0 = par[ini:fin]

            chi += lore(r, D, nu0, self.freq)

        chi += lore(par[-2], par[-1], self.fr, self.freq)

        chi = self.power - chi

        if full:
            return chi
        else:
            return np.sum(chi**2 * self.weigh)

    def make_args(self, method, const, bounds):
        '''
        Returns the proper arguments for the interpolator

        Takes:

        method : the method
        const  : constraints
        biounds: bounds
        '''

        kws = {}

        if method in ['COBYLA', 'SLSQP', 'trust-constr']:
            kws['constraints'] = const

        if method in ['Nelder-Mead', 'Powell',
                      #'COBYLA', 'SLSQP', 'trust-constr',
                      'L-BFGS-B', 'TNC' ]:
            kws['bounds'] =  bounds

        return kws

    def qpo_try_search(self, solution, method='COBYLA'):
        '''
        Fits the QPO searching various frequencies

        It uses the previous solution as initial guesss and adds
        a further Loerntzian with its frequency frozen.
        Everything else is re-fitted.

        It tries a range of frequencies from 1.e-3 up to 10 Hz.
        The returned solution is the one with the minimum chi2.

        Takes:

        solution: the parameters of the solution so far
        method  : the method for the fit

        Returns:

        sola : the new parameters including the QPO
        freco: the frequencies used for the QPO search
        chi2 : list of the chi2 of the frequencies searched
        Qs   : list of the Q values of the frequencies searched
        numQ : the number of the lorentzians after the QPO has been added
        loc  : the index in chi2 and Qs of the solution sola
        '''

        numQ = len(solution) // 3 + 1

        # build the constraints
        # basically:
        # 0 <= r  <= 100
        # 0 <= D  <= 100
        # 0 <= nu <= 100
        # for the QPO:
        # 0 <= r  <= 0.1
        # 0 <= D  <= 1.0
        numo = 3 * numQ - 1
        A  = np.identity(numo)
        lb = np.zeros(numo)
        ub = np.full(numo, 100.)

        lb[-2] = 0.01

        for i in range(numQ - 1):
            ub[3 * i] = 1
        ub[-2] = 0.1
        ub[-1] = 1.0

        #for i in range(1, num):
        #    pos = 3 * (i + 1) - 1
        #    A[pos, pos - 3] = -1

        const = LC(A, lb=lb, ub=ub)

        bounds = np.array([lb, ub]).T

        # the frequencies to try
        freco = np.linspace(1.e-3, 10., 100, endpoint=True)
        #freco = np.geomspace(1.e-3, 10., 100, endpoint=True)


        # these will store the solutions, chi2s and Qs of
        # each trial
        sols = []
        chi2 = []
        Qs   = []
        for self.fr in freco:

            x0 = np.append(solution, [0.085, 1.e-2])
            sols.append(mini(self.funQPO2, x0=x0,
                             **self.make_args(method, const, bounds),
                             method=method, options={'maxiter': 20000}))

            Qs.append(self.fr / sols[-1].x[-1] / 2.)
            chi2.append(self.funQPO2(sols[-1].x))
            print(sols[-1].success)

        # loc of solution
        loc = np.argmin(chi2)

        sola = np.append(sols[loc].x, [freco[loc]])

        return sola, freco, chi2, Qs, numQ, loc

    def qpo_try_max(self, solution, method='COBYLA'):
        '''
        Fits the QPO with a guess from the max of nu f(nu)

        At the first step, the solution is kept frozen, and so
        is the frequency of the QPO. This is the frequency of the max
        of nu f(nu). The free parameters are r and D of the QPO.

        At the second step we combine the solution given, the result
        from the first step and the centroid or peak of the found
        QPO as the initial guess, but refit for everything.

        Takes:

        solution: the parameters of the solution so far
        method  : the method for the fit

        Returns:

        sola : the new parameters including the QPO
        freco: list of the frequency used for the QPO
               (may not be the one in sola) - 1 element
        chi2 : list of the chi2 of the solution - 1 element
        Qs   : list of the Q values of the solution - 1 element
        numQ : the number of the lorentzians after the QPO has been added
        loc  : the index in chi2 and Qs of the solution - 0
        '''

        # the max of nu f(nu)
        loc = np.argmax(self.freq * self.power)

        # for the first step, we need to freeze the
        # given solution and the frequency of the max
        self.res = solution
        self.fr  = self.freq[loc]

        # build the constraints for the QPO:
        # 0 <= r  <= 0.1
        # 0 <= D  <= 100
        A  = np.identity(2)
        lb = np.zeros(2)
        ub = np.full(2, 100.)
        ub[0] = 0.1
        ub[1] = 100.0

        const = LC(A, lb=lb, ub=ub)

        bounds = np.array([lb, ub]).T

        # searches for the best fit at the max of nu f(nu)
        parz = mini(self.funQPO, x0=[0.085, 1.e-2],
                    **self.make_args(method, const, bounds),
                    method=method, options={'maxiter': 20000})
        print(parz.success)

        # calc the guess of nu0 for the second step
        nu0 = np.sqrt(max(self.fr**2 - parz.x[1]**2, 1.e-5))

        # for the second step, everything is fitted,
        # the initial guess is the given solution, the result
        # from the previous step and the nu0 from above
        x0 = np.concatenate((solution, parz.x, [nu0]))

        # refit everything - without constraints
        sola = mini(self.fun, x0=x0, 
                    method=method, options={'maxiter': 20000})

        print(sola.success)
        sola = sola.x

        freco = [self.fr]
        chi2  = [self.fun(sola)]
        Qs    = [sola[-1] / sola[-2] / 2.]
        numQ  = len(sola) // 3
        loc   = 0

        print(numQ)

        return sola, freco, chi2, Qs, numQ, loc

    def fit(self, x0_given=None, plot=False, searchQPO=2,
            ftestlim=0.05, method='COBYLA'):
        '''
        Fits the data, searching for an optimal number of Lorentzians and possibly a QPO

        At least one Lorentzian is fitted, then we try
        an additional one and check with an F-test whether
        the chi2 change is compatible with an improvement.
        If so, we take the extra Lorentzian and check further.

        When we stop, we further add a Lorentzian for the QPO,
        if searchQPO > 0. This search, eventually, refits all
        the Lorentzians.

        Takes:

        x0_given : an initial guess
        plot     : whether to plot the result
        searchQPO: whether to search for an extra Lorentzian
                   that should be the QPO:
                   0: don't
                   1: use qpo_try_search
                   2: use qpo_try_max
        ftestlim : probability to use as stopping condition
                   for the F-test
        method   : the method for the fit

        Returns:

        sola: the parameters of the solution, possibly
              including the QPO if searchQPO > 0
        '''

        if x0_given is None:
            leno = 0
        else:
            leno = len(x0_given) // 3

        num    = 0  # number of Lorentzians for this step
        x0_old = [] # the parameters of previous iteration

        fit_use   = None  # will store the result
        const_use = None  # 
        conv      = False # did we converge?
        cow       = 0     # a counter for the iterations
        chis = []         # stores the chi2 for the F-test
        while not conv:

            # increase the number of Lorentzians
            # and the counter
            num += 1
            cow += 1

            # build constraints
            A  = np.identity(3 * num)
            lb = np.zeros(3 * num)
            ub = np.full(3 * num, 100.)

            for i in range(num):
                ub[3 * i] = 1

            #for i in range(1, num):
            #    pos = 3 * (i + 1) - 1
            #    A[pos, pos - 3] = -1

            const = LC(A, lb=lb, ub=ub)

            bounds = np.array([lb, ub]).T

            # we now build the initial guess for the fitter
            if x0_given is None:

                # if a guess was not given, we just use the
                # previous result plus a deault guess
                x0 = np.append(x0_old, [0.3, 0.5, 0.0])

            else:

                # if a guess was given, we use the old
                # result plus the elements from the given one
                # as long as we can, then we add a default one
                if num <= leno:

                    ini = 3 * (num - 1)
                    x0 = np.append(x0_old, x0_given[ini:ini + 3])

                else:

                    x0 = np.append(x0_old, [0.3, 0.5, 0.0])

            #min_fit = diffevol(self.fun, x0=x0, **self.make_args(method, const, bounds),
            #                   bounds=[[0, 100] for i in range(3 * num)])
            min_fit = mini(self.fun, x0=x0, **self.make_args(method, const, bounds),
                           method=method, options={'maxiter': 20000})

            chis.append(min_fit.fun / (len(self.freq) - 3 * num))

            print(min_fit.success)
            #print('Chi2', chis[-1])

            # now we check for convergence
            # if we did only one step, we cannot check
            if cow == 1:

                conv = False

            else:

                prob = ftest.cdf((chis[-2] - chis[-1]) * (len(self.freq) - 3 * num) / chis[-1] / 3.,
                                 (len(self.freq) - 3 * num + 3),
                                 (len(self.freq) - 3 * num))

                #print(f'Num: {num}, F test: {prob}')

                conv = prob < ftestlim

            if not conv:
                # if we did not converge, store
                # the results so far.
                # convergence means: this step was useless
                fit_use   = min_fit
                const_use = const

            x0_old = min_fit.x

        # we need to bring back the number of Lorentzians by one
        num -= 1


        # if needed, proceed to search for the QPO
        print()
        if   searchQPO == 0:
            sola, numQ = fit_use.x, num
        elif searchQPO == 1:
            sola, freco, chi2, Qs, numQ, loc = self.qpo_try_search(fit_use.x)
        elif searchQPO == 2:
            sola, freco, chi2, Qs, numQ, loc = self.qpo_try_max(fit_use.x)
        print()

        #print(fit_use.x)
        #print(sola)

        # make a plot of the fit showing each lorentzian
        # if asked to do so
        if plot:

            fs = []
            for i in range(numQ):
                ini = 3 * i
                fs.append(lore(*sola[ini:ini + 3], self.freq))

            # if we also searched for QPOs, we show the
            # final residuals and the chi2 and Q value of the search,
            # otherwise we only show the fit
            if searchQPO > 0:
                fig, axi = plt.subplots(3, 1, sharex=True)
            else:
                fig, axi = plt.subplots(1, 1, sharex=True)
                axi = [axi]

            # plot power spectrum with errors
            axi[0].errorbar(self.freq,
                            self.freq * self.power,
                            yerr=self.freq * self.err_power, marker='o')

            # now plot fit

            # total result
            ft = np.copy(fs[0])

            for f in fs[1:]:
                ft += f

            axi[0].plot(self.freq, self.freq * ft, '-')

            # each lorentzian individually
            for f in fs:
                axi[0].plot(self.freq, self.freq * f, ':')

            axi[0].set_xscale('log')
            axi[0].set_yscale('log')
            axi[0].set_ylabel(r'$\nu f(\nu)$')
            axi[-1].set_xlabel(r'$\nu$')

            axi[0].set_ylim([1.e-3, axi[0].get_ylim()[1]])

            if searchQPO > 0:

                # this plots residuals weighted by errors
                axi[1].plot(self.freq, self.fun(sola, full=True) * np.sqrt(self.weigh))
                axi[1].axhline(0, ls='-', color='C1', zorder=-5)
                axi[1].set_xscale('log')
                axi[1].set_ylabel(r'weighted residuals')

                # this plots the chi2 of each frequency we tried
                # and the Q value of the fitted "QPO"
                axi[2].plot(freco, chi2, 'r')
                axi[2].set_ylabel(r'$\chi^2$', color='r')

                ax2 = axi[2].twinx()
                ax2.plot(freco, Qs, 'C0')
                ax2.set_ylabel(r'Q', color='C0')

                for ax in axi:
                    ax.axvline(freco[loc], ls='--', color='k')

            fig.savefig(self.infile + '.pdf', dpi=300)

        return sola


if __name__ == '__main__':

    ext, pre, pos, doQ = '_GX', 'pds_', '', 0   # GX 339-4
    ext, pre, pos,doQ = ''    , '', 'pow', 2   # MAXI J1820 - rebinned
    ext, pre, pos, doQ = '_reb', '', 'pow', 2   # MAXI J1820 - not so rebinned

    # COBYLA   good
    # L-BFGS-B also good
    # SLSQP    not too good I think
    # Powell   not so good either
    meth = 'COBYLA' 

    dira = './maxij1820_powspecs/'

    filelist = sorted(glob.glob(dira + pre + '*' + pos + '.dat'))


    # fit lorentzians to each file
    par_start = None
    pars = []
    for i, infile in enumerate(filelist):

        prima = obs(infile)
        pars.append(prima.fit(par_start, plot=True, searchQPO=doQ,
                              method=meth))


    # write the results to file, one line per data file
    with open('resu' + ext + '.txt', 'w') as filo:

        for i, par in enumerate(pars):

            txt = filelist[i]
            for v in par:
                txt += '  {0: 8.2}'.format(v)
            filo.write(txt + '\n')


    # make a summary figure
    fig, axi = plt.subplots(1, 1)
    axi.set_yscale('log')

    # first let's look for min and max of Q factor
    par  = pars[0]
    Q    = np.array(par[2::3]) / np.array(par[1::3]) / 2.
    mino = Q.min()
    maxo = Q.max()
    
    for i, par in enumerate(pars[1:]):
        Q = abs(np.array(par[2::3]) / np.array(par[1::3]) / 2.)
        mino = min(mino, Q.min())
        maxo = max(maxo, Q.max())

    print('Q lims before limiting', mino, maxo)

    # limit color bar
    mino = max(0., mino)
    maxo = min(2., maxo)

    print('Q lims after  limiting', mino, maxo)

    # now plot, setting colors according to the Q factor
    for i, par in enumerate(pars):

        n = len(par) // 3

        pd = np.array(par[1::3])  # this is the width
        pn = np.array(par[2::3])  # this is nu_0

        cen = np.sqrt(pd**2 + pn**2)

        Qs = abs(pn / pd / 2.)

        sc = axi.scatter(n * [i + 1], cen,
                         marker='o', s=20,
                         c=Qs, cmap='autumn',
                         vmin=mino, vmax=maxo)


    # ok to use the last one, they are all equal
    # in terms of the coloramp scaling used
    cb = plt.colorbar(sc, ax=axi)

    axi.set_xlabel('Obs Num')
    axi.set_ylabel(r'$\nu$')
    cb.set_label('Q', color='r')

    fig.savefig('summa' + ext + '.pdf')

    #input('a')
