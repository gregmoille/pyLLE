import numpy as np
from .analyzedisp import AnalyzeDisp
import scipy.interpolate as itp
import scipy.optimize as optm
import scipy.fftpack as fft
import matplotlib.pyplot as plt
import time
import sys
import subprocess as sub
import os
import inspect
from copy import copy
import pickle as pkl
import shutil
import tempfile
from prettytable import PrettyTable
import warnings
import matplotlib as mpl
import logging
import matplotlib.font_manager as font_manager
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)
    import h5py
import ipdb
path_juliaScript = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/'
tmp_dir = tempfile.gettempdir() + '/'


class LLEsovler(object):
    '''
    Class to solver the Lugiato Lefever Equation
    Initialization input ([]=facultative):

    **res <dict>**

        - Qi <float>: intrinisic Q of the resonator
        - Qc <float>: coupling Q of the resonator
        - R <float>: ring radius
        - gamma <float>: Non linear index of the material

    **sim <dict>**

        - Tscan <float>: length of the simulation (in unit of round trip)
        - mu_fit <list>: number of mode to fit
        - mu_sim <list>: number of mode to simulate
        - dispfile <str> : str pointing to a .mat file where the resonances and corresponding azimuthal mode orders are saved
        - domega_init <float>: initial detuning of the pump 
        - domega_end <float>: final detuning of the pump 
        - [domga_stop] <float>: where to stop the scan in detuning but keep doing the simulation
    
    **debug <bool>**: Save a trace in a logfile in the working directory of the different action pyLLE perform (default = True)
    '''
    _c0 = 299792458
    _ħ = 6.634e-34/(2*np.pi)
    __author__ = "Gregory Moille"
    __copyright__ = "Copyright 2018, NIST"
    __credits__ = ["Gregory Moille",
                    "Qing Li",
                   "Xiyuan Lu",
                   "Kartik Srinivasan"]
    __license__ = "GPL"
    __version__ = "1.0.0"
    __maintainer__ = "Gregory Moille"
    __email__ = "gregory.moille@nist.gov"
    __status__ = "Development"

    def __init__(self, **kwargs):
        self.res = kwargs.get('res', {})
        self.sim = kwargs.get('sim', {})
        self.sim_norm = kwargs.get('sim_norm', None)
        self._debug = kwargs.get('debug', True)

        # find all the needed parameters
        assert 'Qi' in self.res.keys(), 'Please provide Qi'
        assert 'Qc' in self.res.keys(), 'Please provide Qc'
        assert 'R' in self.res.keys(), 'Please provide R'
        assert 'Tscan' in self.sim.keys(), 'Please provide Tscan'
        assert 'dispfile' in self.sim.keys(), 'Please provide dispfile'

        # -- Setup the Logger ---
        if self._debug: 
            FORMATTER = logging.Formatter("[%(asctime)s — LLEsovler.%(funcName)s — %(levelname)s] %(message)s")
            LOG_FILE = "LLE.log"
            open(LOG_FILE, 'a').write('\n' + '-'*75 + '\n')

            self._logger = logging.getLogger(__name__)
            self._logger.setLevel(logging.INFO)
            self._logger.handlers = []
            _loghdl = logging.FileHandler(LOG_FILE)
            _loghdl.setFormatter(FORMATTER)
            self._logger.addHandler(_loghdl)
            self._logger.propagate = False
            self._logger.info('New LLE')
            open(LOG_FILE , 'a').write('-'*75 + '\n')
        else: 
            self._logger = None

    def _Translator(self,D):
        Dnew = {}
        greek ={'α': 'alpha',
                'β':'beta',
                'γ': 'gamma',
                'ω': 'omega',
                'δ': 'd',
                'Δ': 'D',
                'μ': 'mu',
                'λ': 'lambda',
                'ε': 'epsilon',
                'φ': 'phi'}
        for k in D.keys():
            new_k = ''
            for char in list(k):
                if char in list(greek.keys()):
                    new_k += greek[char]
                else:
                    new_k += char
            Dnew[new_k] = D[k]
        return Dnew

    def Analyze(self, plot=False, f=None, ax=None, label=None, plottype='all', zero_lines = True, mu_sim = None):
        '''
        Call pyLLE.analyzedisp.AnalyzeDisp to get the dispersion of the resonator we want to 
        simulate
        '''

        if plot:
            if f is None and ax is None:
                f, ax = plt.subplots()
            elif f is None and not(ax is None):
                if not type(ax) is list:
                    f = ax.figure
                else:
                    f, ax = plt.subplots()
                    print('Only 1 subplots supported, created a new figure')
            elif not(f is None) and ax is None:
                ax = f.axes[0]
            else:
                if type(ax) is list:
                    f, ax = plt.subplots()

        if not('f_pmp' in self.sim.keys()):
            self.sim['f_pmp'] = self._c0/self.sim['lbd_pmp']

        do_plot = plot

        if self._debug:
            Info = '\n\tFilename: {}\n'.format(self.sim['dispfile'])
            Info += '\tf_pmp: {:.3f} THz'.format(self.sim['f_pmp']*1e-12)
            self._logger.info("Analyzing dispersion file" + Info)

        self.sim = self._Translator(self.sim)
        self.res = self._Translator(self.res)

        if mu_sim is None:
            μsim = self.sim['mu_sim']
        else: 
            μsim = mu_sim

        self._analyze = AnalyzeDisp(file=self.sim['dispfile'],
                                    f_center=self.sim['f_pmp'],
                                    rM_fit=self.sim['mu_fit'],
                                    rM_sim=μsim,
                                    R=self.res['R'],
                                    debug=do_plot,
                                    f=f,
                                    ax=ax,
                                    label=label,
                                    plottype=plottype,
                                    zero_lines = zero_lines,
                                    logger = self._logger)
        if do_plot:
            f.canvas.draw()
            plt.pause(0.25)
            f.show()
        PrM_fit, Dint_fit, neff_pmp, ng_pmp = self._analyze.GetDint()
        self._PrM_fit = PrM_fit
        self.sim['Dint'] = Dint_fit
        self.res['ng'] = ng_pmp
        self.res['neff'] = neff_pmp

        self.disp = {}
        self.disp['rf'] = self._analyze.rf
        self.disp['ng'] = self._analyze.ng
        self.disp['neff'] = self._analyze.neff
        self.disp['D'] = self._analyze.D
        self.disp['Dint'] = self._analyze.Dint
        self.disp['dphi'] = self._analyze.dφ
        self.sim['dphi'] = self._analyze.dφ


        self.fDint = f
        self.axDint = ax
        return f, ax 

    def Setup(self):
        '''
        Setup the simulation for the Julia backend. 
        Save the two main dictionary self.sim and self.res into a readable hdf5 file for Julia in the 
        temporary location define by the os
        '''


        # -- Make the hdf5 file --
        # ------------------------------------------------------------
        try:
            os.remove(tmp_dir + 'ParamLLEJulia.h5')
            os.remove(tmp_dir + 'ResultsJulia.h5')
        except:
            pass

        # Check if normalisation or not
        if self.sim_norm == None:
            dic_sim = {'Pin': ('Pin',1e3, 'mW'),
                        'Tscan': ('Tscan',1e-6, 'x1e6 Round Trip'),
                        'f_pmp': ('f_pmp',1e-12, 'THz'),
                        'domega_init': ('δω_init',1e-9/(2*np.pi), 'x2π GHz'),
                        'domega_end': ('δω_end',1e-9/(2*np.pi), 'x2π GHz'),
                        'mu_sim': ('μ_sim',1, ''),
                        'mu_fit': ('μ_fit',1, ''),}
                        
            dic_res = {'R': ('R',1e6, 'µm'),
                        'Qi': ('Qi',1e-6, 'M'),
                        'Qc': ('Qc',1e-6, 'M'),
                        'gamma': ('γ', 1, ''),}
                    

        

            self.sim['debug'] = int(self._debug)
            Info = '-- Solving standard LLE --\n'
            Info += '\tSimulation Parameters\n'
            for k, it in self.res.items():
                if k in dic_res.keys():
                    Info +='\t\t{} = {:.2f} {}\n'.format(dic_res[k][0], it*dic_res[k][1],dic_res[k][2])


            Info += '\tSimulation Parameters\n'
            for k, it in self.sim.items():
                if k in dic_sim.keys():
                    if type(it) is list:
                        Info += '\t\t{} = [{:.2f},{:.2f}] {}\n'.format(dic_sim[k][0], 
                                                                it[0]*dic_sim[k][1],
                                                                it[1]*dic_sim[k][1],
                                                                dic_sim[k][2])
                    else:
                        Info +='\t\t{} = {:.2f} {}\n'.format(dic_sim[k][0], it*dic_sim[k][1],dic_sim[k][2])

            print(Info)
            if self._debug:
                self._logger.info(Info)
                

            # -- create h5file --
            # ipdb.set_trace()
            h5f = h5py.File(tmp_dir + 'ParamLLEJulia.h5', 'w')
            if self._debug:
                self._logger.info('Saving parameters in: {}'.format(tmp_dir + 'ParamLLEJulia.h5'))


            h5f.create_group('sim')
            h5f.create_group('res')
            cnt = 0
            
            for key, it in self.sim.items():    
                if not key == 'δω_disp':
                    if type(it) is str:
                        it = np.string_(it)
                    h5f.create_dataset('sim/{}'.format(key), data=[it])
            for key, it in self.res.items():    
                if not key == 'δω_disp':
                    if type(it) is str:
                        it = np.string_(it)
                    h5f.create_dataset('res/{}'.format(key), data=[it])

            h5f.close()


        else:
            # -- Normalize parameters --
            # ------------------------------------------------------------
            ω0 = 2*np.pi * self.sim['f_pmp']
            ħ = self._ħ
            μ_sim = self.sim['μ_sim']
            c0 = self._c0
            ng = self.res['ng']
            keys = self.sim_norm.keys()
            μ = np.arange(μ_sim[0], μ_sim[-1]+1)
            pmp_ind = np.where(μ==0)[0][0]
            # -- Normalized coupling losses --
            if not('δω_ext' in keys):
                Qc = self.res['Qc']
                Q0 = self.res['Qi']
                self.sim_norm['δω_ext'] = ω0 * (1/Qc)

            # -- Normalized total losses --
            if not('δω_tot' in keys):
                Qc = self.res['Qc']
                Q0 = self.res['Qi']
                self.sim_norm['δω_tot'] = ω0 * \
                    (1/self.res['Qc'] + 1/self.res['Qi'])

            # -- Normalized Pump Power --
            if not('F2' in keys):
                n2 = self.res['n2']
                V = self.res['V']
                Pin = self.sim['Pin']
                δω_ext = self.sim_norm['δω_ext']
                δω_tot = self.sim_norm['δω_tot']
                g0 = n2 * c0 * ħ * ω0**2/(ng**2 * V)
                fact_normE = np.sqrt(8*g0*δω_ext/(δω_tot**3 * ħ * ω0))
                F2 = Pin * fact_normE**2
                self.sim_norm['g0'] = g0
                self.sim_norm['fact_normE'] = fact_normE
                self.sim_norm['F2'] = F2
            
                
            else:
                if not ('fact_normE' in self.sim_norm):
                    n2 = self.res['n2']
                    V = self.res['V']

                    δω_ext = self.sim_norm['δω_ext']
                    δω_tot = self.sim_norm['δω_tot']
                    g0 = n2 * c0 * ħ * ω0**2/(ng**2 * V)
                    fact_normE = np.sqrt(8*g0*δω_ext/(δω_tot**3 * ħ * ω0))
                    self.sim_norm['fact_normE'] = fact_normE
                if not('g0' in keys):
                    n2 = self.res['n2']
                    V = self.res['V']
                    g0 = n2 * c0 * ħ * ω0**2/(ng**2 * V)
                    self.sim_norm['g0'] = g0

            # -- Normalized Detunning --
            if not('α_init' in keys):
                self.sim_norm['α_init'] = -2*self.sim['δω_init']/δω_tot
            if not('α_end' in keys):
                self.sim_norm['α_end'] = -2*self.sim['δω_end']/δω_tot
            if not('α_stop' in keys):
                if not('δω_stop' in self.sim.keys()):
                    self.sim_norm['α_stop'] = "None"
                else:
                    self.sim_norm['α_stop'] = self.sim['δω_stop']
            
            # -- Coupling Dispersion -- 
            if not('δω_disp' in keys):
                if not ('dQc' in self.res.keys()):
                    δω_disp = lambda μ: 1 + 0*μ
                else:
                    δω_disp = lambda μ: 1/self.res.dQc(μ)
            else:
                δω_disp = self.sim_norm['δω_disp']
            
            self.sim_norm['δω_disp'] = δω_disp
            δω_disp_val = δω_disp(μ)
            self.sim_norm['δω_disp_val'] = δω_disp_val
            # -- Misc --
            L =  2*np.pi*self.res['R']
            self.sim_norm['Tscan'] = self.sim['Tscan']
            self.sim_norm['μ_sim'] = self.sim['μ_sim']
            tR = L *self.res['ng']/self._c0
            self.sim_norm['tR'] = tR * self.sim_norm['δω_tot']/2
            self.sim_norm['Dint'] =  self.sim['Dint']
            self.sim_norm['t_end'] = self.sim['Tscan']*self.sim_norm['tR']
            ω = ω0 + μ *2*np.pi/tR
            self.sim_norm['ω'] = ω
            self.sim_norm['μ'] = μ
            # -- Print the Normalized Parameters --
            keys = ['F2','μ_sim',
                  'α_init', 'α_end', 'α_stop',
                  'δω_ext', 'δω_tot',
                  'fact_normE', 'g0',
                  'Tscan','tR','t_end']
            print('– Solving Normalized LLE -')
            print('Normalized Parameters:')
            # for k, it in self.sim_norm.items():
            for k in keys:
                it = self.sim_norm[k]
                if not type(it) is list:
                    try:
                        print('\t{} = {:.2f}'.format(k, it))
                    except:
                        print('\t{} = {}'.format(k, it))
                else:
                    print('\t{} = {}'.format(k, it))

            # -- create h5file --
            # ipdb.set_trace()
            h5f = h5py.File(tmp_dir + 'ParamLLEJulia.h5', 'w')
            h5f.create_group('sim_norm')
            cnt = 0


            
            for key, it in self.sim_norm.items():
                
                if not key == 'δω_disp':
                    if type(it) is str:
                        it = np.string_(it)
                    # print('sim_norm/{}'.format(key))
                    h5f.create_dataset('sim_norm/{}'.format(key), data=[it])
            cnt += 1

            h5f.close()

    def SolveTemporal(self):
        '''
        Call Julia to solver the LLE
        '''
        self._solver = 'temporal'

        if self._debug:
            self._logger.info('Solving Temporal LLE with Julia....')
        print('-'*70)
        date = time.localtime()
        date = "{}-{:0>2}-{:0>2} ".format(date.tm_year,
                                          date.tm_mon,
                                          date.tm_mday) +\
            "{:0>2}:{:0>2}:{:0>2}".format(date.tm_hour,
                                          date.tm_min,
                                          date.tm_sec)
        start = time.time()
        print(date)

        if self.sim_norm == None:
            sub.call(["julia",
                      path_juliaScript + "ComputeLLE.jl", tmp_dir])
        else:
            sub.call(["julia",
                      path_juliaScript + "ComputeLLE_Norm.jl", tmp_dir])
        time_taken = time.time() - start
        end = time.time()
        hours, rem = divmod(end-start, 3600)
        minutes, seconds = divmod(rem, 60)
        time_taken = "Simulation Time " + \
            "{:0>2}h:{:0>2}min:{:0>4.1f}s".format(int(hours),
                                                  int(minutes),
                                                  seconds)
        print('\n')
        print(time_taken)
        print('-'*70)

    def SolveSteadySteate(self):
        '''
        Newton Method to find the root of the steady state equation
        '''

        self._solver = 'steady'
        if self._debug:
            self._logger.info('Solving Steady State LLE with Python....')
        print('-'*70)
        date = time.localtime()
        date = "{}-{:0>2}-{:0>2} ".format(date.tm_year,
                                          date.tm_mon,
                                          date.tm_mday) +\
            "{:0>2}:{:0>2}:{:0>2}".format(date.tm_hour,
                                          date.tm_min,
                                          date.tm_sec)
        start = time.time()
        print(date)

        # -- CHeck Parity of the µ --
        μ_sim = self.sim['mu_sim']
        ind = 0
        if not μ_sim[0] == -μ_sim[1]:
            μ_sim = [-np.max(np.abs(μ_sim)), np.max(np.abs(μ_sim))]
            INFO = 'Not symmetric mode calculation -> switching to it with µ_sim = {}\n'.format(μ_sim)
            print(INFO)
            if self._debug:
                self._logger.info(INFO)
            self.Analyze(mu_sim = μ_sim)


        # -- setting up parameters -- 
        β2 = self._analyze.β2
        Pin = self.sim['Pin']
        γ = self.res['gamma']
        L = 2*np.pi*self.res['R']
        ω0 = self.sim['f_pmp']*2*np.pi
        Q0 = self.res['Qi']
        Qc = self.res['Qc']
        tR = L*self.res['ng']/self._c0
        α = 1/2 * (ω0/Q0 + ω0/Qc) * tR
        θ = ω0/Qc*tR 
        δω = -self.sim['domega']* tR
        
        nlc = -1j*γ*L
        μ = np.arange(μ_sim[0], μ_sim[1]+1)
        pmp_ind = np.where(μ == 0)[0][0]
        ω = μ*2*np.pi/tR + ω0
        ν = ω/(2*np.pi)
        # -- setting up input power -- 
        Ein = np.zeros(μ.size)
        Ein[pmp_ind] = np.sqrt(Pin) * μ.size
        Ein_couple = np.sqrt(θ)*Ein
    
        # -- Define Initial Guess --
        sech = lambda x: 1/np.cosh(x)
        η = δω/α # need to be fixed
        B0 = np.sqrt(2*η)
        f = np.sqrt(θ* Pin* γ* L/α**3)
        φ0 = np.arccos(np.sqrt(8*η)/np.pi/f)
        τ = np.linspace(-0.5,0.5, μ.size)*tR
        Φ0 =  f/η**2 -1j* f/η
        ut0 = np.sqrt(α/γ/L) * (Φ0+ B0 * np.exp(1j*φ0) * sech(B0*np.sqrt(α/(np.abs(β2)*L))*τ))
        Em0 = fft.fftshift(fft.fft(ut0))
        x_init = np.concatenate((Em0.real, -Em0.imag))
        
        # -- Define the Steady State LLE Equation --
        φ = -α + 1j*δω - 1j*self.sim["dphi"]
        Em= lambda xx:  xx[0:int(xx.shape[0]/2)] + 1j*xx[int(xx.shape[0]/2)]
        Ut=  lambda xx:  fft.ifft(Em(xx));
        fm= lambda xx: φ*Em(xx) + nlc*fft.fft(np.abs(Ut(xx))**2*Ut(xx)) + Ein_couple;
        fvec= lambda xx: np.concatenate((fm(xx).real, fm(xx).imag))

        # -- Solver the Steady State -- 
        out = optm.root(fvec, x_init, method='lm', jac=None, tol=1e-20)
        Ering = Em(out.x)/μ.size
        Ewg = Ein/μ.size -Ering*np.sqrt(θ)


        
        f, ax = plt.subplots()
        ax.plot(1e-12*ν, 20*np.log10(np.abs(Ering)), label='Ring')
        ax.plot(1e-12*ν, 20*np.log10(np.abs(Ewg)), label='Waveguide')
        ax.legend()
        f.show()

        return Ering, Ewg, f, ax

    def RetrieveData(self):
        '''
        Load the output hdf5 save by julia and transform it in a user-friendly dictionary to be more pyhtonistic 
        '''


        drct = tmp_dir
        S = h5py.File(tmp_dir + 'ResultsJulia.h5', 'r')
        if self._debug:
            self._logger.info('Retrieving results from Julia in {}'.format(tmp_dir + 'ResultsJulia.h5'))
        sol = {}
        keys = ['u_probe','Em_probe', 'Ewg']
        for k in keys:
            rl = 'Results/{}Real'.format(k)
            im = 'Results/{}Imag'.format(k)
            sol[k] = S[rl][:] + 1j*S[im][:]
        sol['ω'] = S['Results/ωReal'][:]
        sol['comb_power'] = S['Results/comb_powerReal'][:]
        sol['detuning'] = S["Results/detuningReal"][:]
        S.close()
        os.remove(tmp_dir + 'ParamLLEJulia.h5')
        os.remove(tmp_dir + 'ResultsJulia.h5')
        sol['freq'] = sol['ω']/(2*np.pi)

        # Normalize to SI
        # g0 = self.sim_norm['g0'] 
        # δω_tot = self.sim_norm['δω_tot'] 
        # # normE = np.sqrt(2*g0/δω_tot)
        # normE = self.sim_norm['fact_normE']
        # beta = self.sim_norm['Dint']/δω_tot
        # Nprobe = sol['Em_probe'].shape[1]
        # t_end = self.sim_norm['t_end']
        # μ_sim = self.sim['μ_sim']

        # tau = np.linspace(0,t_end, Nprobe)

        # μ = np.arange(μ_sim[0], μ_sim[-1]+1)


        # cnt = 0
        # for tt in tau: 
        #     phi = np.exp(-1j*beta*tt)
        #     Epb = sol['Em_probe'][:,cnt]
        #     Ewg = sol['Ewg'][:,cnt]
        #     sol['Em_probe'][:,cnt] = phi*np.conj(Epb)/normE
        #     sol['Ewg'][:,cnt] = phi*np.conj(Ewg)/normE
        #     th_cnt = 0
        sol['theta'] = np.linspace(-np.pi,np.pi, sol['u_probe'].shape[0])
        # sol['u_probe'] = np.conj(sol['u_probe'])
        self.sol = sol

    def PlotCombPower(self):
        '''
        Plot a figure with 3 subplots. 

        - Top subplot = map of the spectra for the steps taken by the LLE (step sub-sampled to be 1000)
        - middle subplot = temporal map of the intensity inside the resonator for the steps of the LLE
        - bottom subplot = normalized comb power

        **Output**

            - f, ax handle of figure and axes of the matplotlib figure displayed
        '''
        cmap = plt.get_cmap("tab10")
        freq = self.sol['freq']*1e-12
        Epb = self.sol['Ewg']
        Epb[Epb==0] = 1e-20
        Epb = 10*np.log10(np.abs(Epb)**2)
        Epb = Epb - Epb.max()
        # Epb[Epb<=0] = 1e-100
        f, ax = plt.subplots(3, 1, sharex=True)
        ax[0].pcolormesh(np.arange(0, 1000), freq,Epb,
                         rasterized=True, 
                         vmin = -100)
        theta = self.sol['theta']
        ax[1].pcolormesh(np.arange(0, 1000), theta/np.pi,
                         (np.abs(self.sol['u_probe'])**2),
                         rasterized=True, vmin=0)
        ax[2].plot(np.arange(0, 1000), self.sol['comb_power'] /
                   self.sol['comb_power'].max())
        ax = np.append(ax, ax[2].twinx())
        ax[3].plot(np.arange(0, 1000),self.sol['detuning']*1e-9/(2*np.pi),
                    c = cmap.colors[1])
        ax[0].set_ylabel('Frequency (THz)')
        ax[1].set_ylabel('Angle (x π)')
        ax[2].set_xlabel('LLE Step (sub-sampled)')
        ax[2].set_ylabel('Norm. Comb Pwr')
        ax[3].set_ylabel('Detuning (GHz)')
        ax[0].set_xlim([0,1000])
        f.show()

        self.fPcomb = f
        self.axPcomb = ax

        return f, ax

    def PlotCombSpectra(self, ind, f=None, ax=None, label=None, pwr='both'):
        '''
        Plot the spectra for a given index in the 1000 sub-sampled LLE step 

        **Input** 

            - ind <ind>: index in the LLE step to plot the spectra
            - f <obj>:  matplotlib figure handle (if None, new figure)
            - ax <obj>: matplotlib axe handle 
            - label <str>: label for the legend
            - pwr <str>: 'both', 'ring', 'wg' depending on the spectra wanted (inside the ring, the waveguide or both)

        **Output**

            - f <obj>:  matplotlib figure handle
            - ax <obj>: matplotlib axe handle
            - freq <numpy.array>: frequency in Hz
            - Sout <numpy.array>: spectral density of power in the waveguide (dBm)
            - Sring <numpy.array>: spectral density of power in the ring (dBm)

        '''

        freq = self.sol['freq']*1e-12
        if f is None and ax is None:
            f, ax = plt.subplots()
        elif f is None and not(ax is None):
            if not type(ax) is list:
                f = ax.figure
            else:
                f, ax = plt.subplots()
                print('Only 1 subplots supported, created a new figure')
        elif not(f is None) and ax is None:
            ax = f.axes[0]
        else:
            if type(ax) is list:
                f, ax = plt.subplots()

        Sring = 10*np.log10(np.abs(self.sol['Em_probe'][:, ind])**2)
        Sout = 10*np.log10(np.abs(self.sol['Ewg'][:, ind])**2)

        if pwr.lower() == 'both':
            ax.plot(freq, Sout, label='Output P')
            ax.plot(freq, Sring, '.', ms=4, label='In ring P')
        if pwr.lower() == 'ring':
            ax.plot(freq, Sring, ms=4, label=label)
        if pwr.lower() == 'wg':
            ax.plot(freq, Sout, label=label)

        ax.set_ylabel('Power (dBm)')
        ax.set_xlabel('Frequency (THz)')
        if not(label is None):
            ax.legend()

        f.canvas.draw()
        f.show()
        plt.pause(0.25)

        self.fSpectra = f
        self.axSpectra = ax

        return freq, Sout, Sring, f, ax

    def PlotSolitonTime(self, ind):
        freq = self.sol['freq']
        tR = 2*np.pi*self.res['R']*self.res['ng']/self._c0
        f, ax = plt.subplots()
        τ = np.linspace(-0.5, 0.5, freq.size) * tR
        U = np.abs(self.sol['u_probe'][:,ind])**2
        ax.plot(τ*1e12 , U/U.max())
        ax.set_xlabel('Time (ps)')
        ax.set_ylabel('Soliton Energy (a.u)')
        f.show()

        return τ, U, f, ax

    def SaveResults(self, fname, path='./'):
        '''
        Save the whole class with pickle to be able to easilly call it back or retrieve the results after saving

        **Input**

            - fname <str>: name to save. The '.pkl' extension will be added
            - path <str>: path to save the results (defaults './')
        '''


        to_save = copy(self)
        to_save.sim_norm.pop('domega_disp', None)
        to_save.sim.pop('domega_disp', None)
        fname = path + fname + '.pkl'
        print(fname)
        pkl.dump(to_save, open(fname,'bw'))

    def __repr__(self):
        to_print = ''
        to_print = 'Dispersion load from:\t{}\n\n'.format(self.sim['dispfile']) 
        to_print += 'Resonator Parameters:\n'
        res_table = PrettyTable(['Parameters', 'Value', 'Units'])
        res_table.add_row(['R', "{:.3f}".format(self.res['R']*1e6),'µm'])
        res_table.add_row(['Qi', "{:.3f}".format(self.res['Qi']*1e-6),'x1e6'])
        res_table.add_row(['Qc', "{:.3f}".format(self.res['Qc']*1e-6),'x1e6'])
        if 'gamma' in self.res:
            res_table.add_row(['γ', "{:.3f}".format(self.res['gamma']),''])
        if 'n2' in self.res:
            res_table.add_row(['n2', "{:.3f}".format(self.res['n2']*1e19),'x1e-19 m2/W'])
        to_print += res_table.get_string()
        to_print += '\n'

        to_print += 'Simulation Parameters:\n'
        sim_table = PrettyTable(['Parameters', 'Value', 'Units'])

        if 'Pin' in self.sim:
            sim_table.add_row(['Pin',"{:.3f}".format(self.sim['Pin']*1e3),'mW'])
        if 'f_pmp' in self.sim:
            sim_table.add_row(['f_pmp',"{:.3f}".format(self.sim['f_pmp']*1e-12),'THz'])
        if 'μ_sim' in self.sim:
            sim_table.add_row(['μ_sim',"{}".format(self.sim['mu_sim']),''])
        if 'Tscan' in self.sim:
            sim_table.add_row(['Tscan',"{:.3f}".format(self.sim['Tscan']*1e-5),'x1e5 tR'])        
        if 'domega_init' in self.sim:
            sim_table.add_row(['δω_init',"{:.3f}".format(self.sim['domega_init']*1e-9),'GHz'])
        if 'domega_end' in self.sim:
            sim_table.add_row(['δω_end',"{:.3f}".format(self.sim['domega_end']*1e-9),'GHz'])
        to_print += sim_table.get_string()
        to_print += '\n'

        # to_print += 'Normalized Simulation Parameters:\n'
        # aa = self.sim_norm['δω_disp']
        # fun = str(inspect.getsourcelines(aa)[0]).split('lambda')[-1].split(':')[-1].strip('["\\n"]')[1::]
        # sim_norm_table = PrettyTable(['Parameters', 'Value'])
        # sim_norm_table.add_row(['F2',"{:.3f}".format(self.sim_norm['F2'])])
        # sim_norm_table.add_row(['α_init',"{:.3f}".format(self.sim_norm['α_init'])])
        # sim_norm_table.add_row(['α_end',"{:.3f}".format(self.sim_norm['α_end'])])
        # sim_norm_table.add_row(['α_stop',"{}".format(self.sim_norm['α_stop'])])
        # sim_norm_table.add_row(['δω_disp',"{}".format(fun)])
        # to_print += sim_norm_table.get_string()
        # to_print += '\n'
        
        return to_print
