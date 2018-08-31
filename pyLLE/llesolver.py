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
import matplotlib.gridspec as gridspec
from matplotlib import ticker
import matplotlib.font_manager as font_manager
from datetime import datetime
from plotly import tools
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
import plotly.graph_objs as go

with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)
    import h5py



backend = mpl.get_backend()
path_juliaScript = os.path.dirname(os.path.abspath(__file__))
path_juliaScript = os.path.join(path_juliaScript, 'ComputeLLE.jl')
tmp_dir = tempfile.mkdtemp()
os.rmdir(tmp_dir) # delete folder as the temp file will not be in it but will be used as a prefix
print('-'*50)
print('Path to temporary dir to save .h5 files with prefix:')
print(tmp_dir)
print('-'*50)
print('Path to Julia script: ')
print(path_juliaScript)
print('-'*50)



# Check which type of python we are launching

try:
    className = get_ipython().__class__.__name__
    if className == 'ZMQInteractiveShell':
        pyType = 'jupyter'
    elif className == 'TerminalInteractiveShell':
        pyType = 'ipython' 
except:
    # launching trhough a normal python
    pyType = 'normal'

# print(pyType)


class MyLogger():
    '''
    Custom made logger as the logger default package cannot be pickled
    '''

    def __init__(self, fname):
        self.fname = fname
        open(self.fname,'a').write('\n' + '-'*75 + '\n')
    
    def info(self, method, message):
        time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        mess = '[ ' + time + ' - ' + method + ' ] ' + message + '\n'
        open(self.fname,'a').write(mess)
    
class Latexify():
    '''
    Class that handle saving the figures in a nice way compatible with 
    the column/page size of different latex template

    input [] = optional: 
        - figname = name to save the figure (without extension)
        -
        - fig = matplotlib handle to the figure
        - [fig_width]: default = 1column
        - [frmt]: default = pdf
        - [fig_height] : default = 6.5
        - [font_size] : default = 8
    '''
    __author__ = "Gregory Moille"
    __copyright__ = "Copyright 2018, NIST"
    __credits__ = ["Gregory Moille",
                   "Kartik Srinivasan"]
    __license__ = "GPL"
    __version__ = "2.0.1"
    __maintainer__ = "Gregory Moille"
    __email__ = "gregory.moille@nist.gov"
    __status__ = "Development"

    def __init__(self, **kwarg):
        # get parameters
        figname = kwarg.get('figname', '')
        fig = kwarg.get('fig', None)
        fig_width = kwarg.get('fig_width', '1column')
        self.frmt = kwarg.get('frmt', 'pdf')
        self.fig_height = kwarg.get('fig_height', 6.5)
        self.font_size = kwarg.get('font_size', 8)
        if isinstance(fig_width, str):
            if fig_width.lower() == '1column':
                self.fig_width = 8.6
            elif fig_width.lower() == '2column':
                self.fig_width = 14
            if fig_width.lower() == '1columnbeamer':
                self.fig_width = 10.79846 / 2
                if not kwarg.get('fig_height', False):
                    self.fig_height = 6.5 * 10.79846 / 24
            if fig_width.lower() == '2columnbeamer':
                self.fig_width = 10.79846
                if not kwarg.get('fig_height', False):
                    self.fig_height = 6.5 * 10.79846 / 14
        else:
            self.fig_width = fig_width

        inch = 2.54
        self.fig_width = self.fig_width / inch
        self.fig_height = self.fig_height / inch
        self.f = fig
        self.figname = figname
        self.SavePlot()

    def SavePlot(self):
        # -- Define Font Properties --
        # -----------------------------------------------------
        # fontpath = '/System/Library/Fonts'
        font_prop = font_manager.FontProperties(size=8)
        plt.ion()
        # -- Define Param of Plot --
        # -----------------------------------------------------
        params = {'backend': 'ps',
                  'text.latex.preamble': [r'\usepackage{gensymb}',
                                          r'\usepackage{siunitx}',
                                          r'\sisetup{detect-all}',
                                          r'\usepackage{helvet}',
                                          r'\usepackage{sansmath}',
                                          r'\sansmath', ],
                  'text.latex.unicode': False,
                    # fontsize for x and y labels(was 10)
                  'axes.labelsize': self.font_size,
                  'axes.titlesize': self.font_size,
                  'axes.linewidth': 0.5,
                  'xtick.major.width': 1,
                  'xtick.minor.width': 1,
                  'ytick.major.width': 1,
                  'ytick.minor.width': 1,
                  'legend.fontsize': self.font_size,  # was 10
                  'xtick.labelsize': self.font_size,
                  'ytick.labelsize': self.font_size,
                  'text.usetex': False,
                  'figure.figsize': [self.fig_width, self.fig_height],
                  # 'font.family': 'sans-serif',
                  'font.size': self.font_size,
                  'lines.linewidth': 0.25,
                  }
        mpl.rcParams.update(params)

        plt.pause(0.1)
        self.f.set_facecolor('None')
        # -- Redo Font --
        # -----------------------------------------------------
        for ax in self.f.axes:
            ax.xaxis.set_ticks_position('both')
            ax.yaxis.set_ticks_position('both')

            for line in ax.yaxis.get_ticklines():
                line.set_markeredgewidth(0.25)
                line.set_markersize(3)
            plt.tick_params(which='minor', length=2, width=0.25)

            for line in ax.xaxis.get_ticklines():
                line.set_markeredgewidth(0.25)
                line.set_markersize(3)
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(self.font_size) 
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(self.font_size) 
            for axis in ['top', 'bottom', 'left', 'right']:
                ax.spines[axis].set_visible(True)
                ax.spines[axis].set_edgecolor('k')
                ax.spines[axis].set_linewidth(0.25)
            ax.axesPatch.set_facecolor('None')  
            
            # ax.grid('off')
            xstr = ax.get_xlabel()
            ax.set_xlabel(xstr, size=self.font_size)
            ystr = ax.get_ylabel()
            ax.set_ylabel(ystr, size=self.font_size)
            leg = ax.axes.get_legend()
            # Check if ther is a legend
            # ----------------------------------------
            if leg:
                txt = leg.properties()['texts']
                lbl_lines = leg.properties()['lines']
                for ii in txt:
                    ii.set_fontsize(self.font_size)
                for ii in lbl_lines:
                    ii.set_linewidth(0.25)
            #change linestyle for line2D
            for ch in ax.get_children():
              if type(ch) == mpl.lines.Line2D:
                ch.set_linewidth(0.25)
                ch.set_markersize(2)
                ch.set_markeredgewidth(4)
                ch.set_markeredgecolor('None')
        # -- Update plot --
        # -----------------------------------------------------
        self.f.set_size_inches(self.fig_width, self.fig_height)
        self.f.canvas.draw()
        plt.draw()
        plt.pause(0.05)
        if not self.figname == '':
            self.f.savefig(self.figname + '.' + self.frmt, format = self.frmt)
        plt.pause(0.05)
        self.GoBackToNormal()
        plt.ioff()

    def GoBackToNormal(self):
        file = mpl.get_data_path() + '/matplotlibrc'
        mpl.rcParams.update(mpl.rc_params_from_file(
            file))
        self.f.canvas.draw()

class LLEsovler(object):
    '''
    Class to solve the Lugiato Lefever Equation
    Initialization input ([]=facultative):

    **res <dict>**

        - Qi <float>: intrinsic Q of the resonator
        - Qc <float>: coupling Q of the resonator
        - R <float>: ring radius
        - gamma <float>: Non linear index of the material

    **sim <dict>**

        - Tscan <float>: length of the simulation (in unit of round trip)
        - mu_fit <list>: number of mode to fit
        - mu_sim <list>: number of mode to simulate
        - dispfile <str> : str pointing to a .csv file where the azimuthal mode orders and corresponding resonances are saved
        - domega_init <float>: initial detuning of the pump 
        - domega_end <float>: final detuning of the pump 
        - [domga_stop] <float>: where to stop the scan in detuning but keep doing the simulation
    
    **debug <bool>**: Save a trace in a logfile in the working directory of the different actions pyLLE perform (default = True)
    '''
    _c0 = 299792458
    hbar = 6.634e-34/(2*np.pi)
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
        self._plotPower = None
        self._plotSpecta = None
        self._indSpectra = 0
        self._plotTime = None
        self._indTime = 0

        # find all the needed parameters
        assert 'Qi' in self.res.keys(), 'Please provide Qi'
        assert 'Qc' in self.res.keys(), 'Please provide Qc'
        assert 'R' in self.res.keys(), 'Please provide R'
        assert 'Tscan' in self.sim.keys(), 'Please provide Tscan'
        assert 'dispfile' in self.sim.keys(), 'Please provide dispfile'

        # -- Setup the Logger ---
        if self._debug: 
            self._logger = MyLogger("LLE.log")
            self._logger.info('__init__', 'New LLE')
        else: 
            self._logger = None

    def _Translator(self,D):
        Dnew = {}
        self._greek ={'α': 'alpha',
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
                if char in list(self._greek.keys()):
                    new_k += self._greek[char]
                else:
                    new_k += char
            Dnew[new_k] = D[k]
        return Dnew

    def Analyze(self, plot=False, f=None, ax=None, label=None, plottype='all', zero_lines = True, mu_sim = None):
        '''
        Call pyLLE.analyzedisp.AnalyzeDisp to get the dispersion of the resonator we want to simulate
        '''

        if plot and (not pyType == 'jupyter'):
            if f is None and ax is None:
                f, ax = plt.subplots(dpi=120)
            elif f is None and not(ax is None):
                if not type(ax) is list:
                    f = ax.figure
                else:
                    f, ax = plt.subplots(dpi=120)
                    print('Only 1 subplots supported, created a new figure')
            elif not(f is None) and ax is None:
                ax = f.axes[0]
            else:
                if type(ax) is list:
                    f, ax = plt.subplots(dpi=120)
        else:
            f = None
            ax = None
        if not('f_pmp' in self.sim.keys()):
            self.sim['f_pmp'] = self._c0/self.sim['lbd_pmp']

        do_plot = plot

        if self._debug:
            Info = '\n\tFilename: {}\n'.format(self.sim['dispfile'])
            Info += '\tf_pmp: {:.3f} THz'.format(self.sim['f_pmp']*1e-12)
            self._logger.info('LLEsovler.Analyze', "Analyzing dispersion file" + Info)

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
                                    logger = self._logger,
                                    pyType = pyType)
        if do_plot and (not pyType == 'jupyter'):
            f.canvas.draw()
            plt.pause(0.25)
            f.show()

        PrM_fit, Dint_fit, neff_pmp, ng_pmp, f, ax = self._analyze.GetDint()
        self._PrM_fit = PrM_fit
        self.sim['Dint'] = Dint_fit
        self.res['ng'] = ng_pmp
        self.res['neff'] = neff_pmp

        self.disp = {}
        self.disp['freq'] = self._analyze.rf
        self.disp['ng'] = self._analyze.ng
        self.disp['neff'] = self._analyze.neff
        self.disp['D'] = self._analyze.D
        self.disp['Dint'] = self._analyze.Dint
        self.disp['dphi'] = self._analyze.dφ
        self.sim['dphi'] = self._analyze.dφ

        if (not pyType == 'jupyter'):
            self.fDint = f
            self.axDint = ax
            return f, ax 

    def Setup(self):
        '''
        Setup the simulation for the Julia back-end. 
        Save the two main dictionary self.sim and self.res into a readable hdf5 file for Julia in the temporary location define by the os
        '''
        # -- Make the hdf5 file --
        # ------------------------------------------------------------
        try:
            os.remove(tmp_dir + 'ParamLLEJulia.h5')
            os.remove(tmp_dir + 'ResultsJulia.h5')
        except:
            pass

        dic_sim = {'Pin': ('Pin',1e3, 'mW'),
                    'Tscan': ('Tscan',1e-6, 'x1e6 Round Trip'),
                    'f_pmp': ('f_pmp',1e-12, 'THz'),
                    'domega_init': (u'\u03B4\u03C9_init',1e-9/(2*np.pi), u'x2\u03C0 GHz'),
                    'domega_end': (u'\u03B4\u03C9_end',1e-9/(2*np.pi), u'x2\u03C0 GHz'),
                    'mu_sim': (u'\u03BC_sim',1, ''),
                    'mu_fit': (u'\u03BC_fit',1, ''),}
                    
        dic_res = {'R': ('R',1e6, 'µm'),
                    'Qi': ('Qi',1e-6, 'M'),
                    'Qc': ('Qc',1e-6, 'M'),
                    'gamma': (u'\u03B3', 1, ''),}
                

    

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
            try:
                self._logger.info('LLEsovler.Setup', Info)
            except: 
                Info = ''.join([self._greek[ii] if ii in self._greek.keys() else ii for ii in Info])
                self._logger.info('LLEsovler.Setup', Info)
            

        # -- create h5file --
        # ipdb.set_trace()
        h5f = h5py.File(tmp_dir + 'ParamLLEJulia.h5', 'w')
        if self._debug:
            self._logger.info('LLEsovler.Setup','Saving parameters in: {}'.format(tmp_dir + 'ParamLLEJulia.h5'))


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

    def SolveTemporal(self, tol = 1e-3, maxiter = 6, step_factor = 0.1):
        '''
        Call Julia to solve the LLE
        '''
        self._solver = 'temporal'

        if self._debug:
            self._logger.info('LLEsovler.SolveTemporal','Solving Temporal LLE with Julia....')
            Info = 'tol = {} -- maxiter = {} step_factpr = {}'.format(tol, maxiter, step_factor)
            self._logger.info('LLEsovler.SolveTemporal',Info)


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


        if sys.platform == 'darwin':
            julia = 'julia'
        if sys.platform == 'linux2':
            julia = 'julia'
        if sys.platform == 'win32':
            julia = os.path.expanduser('~') + '\\AppData\\Local\\Julia-0.6.4\\bin\\julia.exe'

        # if self.sim_norm == None:
        command = [julia, path_juliaScript , tmp_dir, str(tol), str(maxiter), str(step_factor)]
        self.JuliaSolver = sub.Popen(command, stdout=sub.PIPE)
        print('Launching Julia')
        line = ''
        len_lin = len(line)
        fname = tmp_dir + 'log.log' 
        conv_err = False

        def Pbar(perc, pgrs, tb_up, prev_line):
            bar_width = 50
            pgrs = '•'*int(np.floor(perc/2))
            width = ' '* (bar_width - len(pgrs))
            perc_str = ' {}%'.format(perc)
            line = 'Computing LLE [' + pgrs  + width  + ']' +  perc_str
            if conv_err:
                line = line + ' /!\ Convergence issue'
                if self._debug:
                    self._logger.info('LLEsovler.SolveTemporal','/!\ Convergence issue')
            length = len(line)
            return line, length, pgrs, tb_up


        tb_up = 2
        pgrs = ''
        prev_line = ''
        perc_old = 0
        perc = -1
        
        line = ''

        while not perc == 100 and self.JuliaSolver.poll() == None:
            try:
                ll = open(fname).readlines()[-1].strip()
                try: 
                    perc = int(ll)  
                    if not perc_old == perc:
                        line, length, pgrs, tb_up= Pbar(perc, pgrs, tb_up, prev_line)
                        prev_line = line
                        print(line, end="\r")
                        perc_old = perc

                except Exception as e:
                    if ll.split()[0] == 'Failed':
                        conv_err = True
            except:
                pass
        
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
        if self._debug:
            self._logger.info('LLEsovler.SolveTemporal', time_taken)
            

    def SolveSteadySteate(self):
        '''
        Newton Method to find the root of the steady state equation
        '''

        self._solver = 'steady'
        if self._debug:
            self._logger.info('LLEsovler.SolveSteadySteate','Solving Steady State LLE with Python....')
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
        μ_sim = copy(self.sim['mu_sim'])
        ind = 0
        if not μ_sim[0] == -μ_sim[1]:
            μ_sim = [-np.max(np.abs(μ_sim)), np.max(np.abs(μ_sim))]
            INFO = 'Not symmetric mode calculation -> switching to it with µ_sim = {}\n'.format(μ_sim)
            print(INFO)
            if self._debug:
                try:
                    self._logger.info('LLEsovler.SolveSteadySteate',INFO)
                except:
                    Info = ''.join([self._greek[ii] if ii in self._greek.keys() else ii for ii in Info])
                    self._logger.info('LLEsovler.SolveSteadySteate',INFO)

            dum = self.Analyze(mu_sim = μ_sim, plot = False, f = None)


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

        self.steady = {'Ering':Ering, 'Ewg':Ewg}
        if not pyType == 'jupyter':
            f, ax = plt.subplots(dpi=120)
            ax.plot(1e-12*ν, 20*np.log10(np.abs(Ering)), label='Ring')
            ax.plot(1e-12*ν, 20*np.log10(np.abs(Ewg)), label='Waveguide')
            ax.legend()
            f.show()
            return f, ax
        else:
            trace0 = go.Scatter(x = 1e-12*ν,y = 20*np.log10(np.abs(Ering)),
                            mode = 'lines', name='Res. Power')
            trace1 = go.Scatter(x = 1e-12*ν,y = 20*np.log10(np.abs(Ewg)),
                            mode = 'lines', name='Out Power')
            data = [trace0, trace1]
            layout = dict(xaxis = dict(title = 'Frequency (THz)'),
                  yaxis = dict(title = 'Power (dBm)'),
                  )
            fig = go.Figure(data=data, layout=layout)
            iplot(fig)

    def RetrieveData(self):
        '''
        Load the output hdf5 saved by julia and transform it in a user-friendly dictionary to be more pythonistic 
        '''

        time.sleep(0.5)
        drct = tmp_dir
        S = h5py.File(tmp_dir + 'ResultsJulia.h5', 'r')
        if self._debug:
            self._logger.info('LLEsovler.RetrieveData','Retrieving results from Julia in {}'.format(tmp_dir + 'ResultsJulia.h5'))
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
        sol['theta'] = np.linspace(-np.pi,np.pi, sol['u_probe'].shape[0]) 
        self.sol = sol

    def PlotCombPower(self, do_matplotlib = False):
        '''
        Plot a figure with 3 subplots. 

        - Top subplot = map of the spectra for the steps taken by the LLE (step sub-sampled to be 1000)
        - middle subplot = temporal map of the intensity inside the resonator for the steps of the LLE
        - bottom subplot = normalized comb power

        **Output**

            - f, ax:  handle of figure and axes of the matplotlib figure displayed
        '''


        freq = self.sol['freq']*1e-12
        Epb = self.sol['Ewg']
        Epb[Epb==0] = 1e-20
        Epb = 30 + 10*np.log10(np.abs(Epb)**2)

        E2 = (np.abs(self.sol['u_probe'])**2)
        E2 = E2/E2.max()
        tR = 2*np.pi*self.res['R']*self.res['ng']/self._c0
        t = np.linspace(-0.5, 0.5, freq.size) * tR

        CmbPow = self.sol['comb_power'] /self.sol['comb_power'].max()
        det = self.sol['detuning']*1e-9/(2*np.pi)

        step = np.arange(0, 1000)

        if not pyType == 'jupyter' or do_matplotlib:
            # --  Create the Figure -- 
            f = plt.figure()
            gs = gridspec.GridSpec(3,2, width_ratios=[1,0.015],wspace=0.05)
            ax = [None]*6
            ax[0] = plt.subplot(gs[0])
            ax[1] = plt.subplot(gs[1])
            ax[2] = plt.subplot(gs[2],sharex=ax[0])
            ax[3] = plt.subplot(gs[3])
            ax[4] = plt.subplot(gs[4],sharex=ax[0])
            cmap = plt.get_cmap("tab10")

            # -- Plot Everything -- 
            aa = ax[0].pcolormesh(step, freq,Epb,
                             rasterized=True, 
                             vmin = Epb.max()-120,
                             vmax = Epb.max())

            bb = ax[2].imshow(E2,aspect='auto',
                    origin = 'lower', 
                    interpolation='bessel')
            tr_12 = 1e-12*np.floor(1e12*tR)/2
            # print(np.argmin(np.abs(tr_12- t)))
            ind = [np.argmin(np.abs(-tr_12- t)),
                   np.argmin(np.abs(t)),
                   np.argmin(np.abs(tr_12- t))]
            ax[2].set_yticks(ind)
            ax[2].set_yticklabels([-tr_12*1e12,
                                     0,
                                    tr_12*1e12])


            ax[4].plot(step, CmbPow)
            ax.append(ax[4].twinx())
            ax[6].plot(step,det,
                        c = cmap.colors[1])

            # -- Make it Pretty --
            ax[0].set_ylabel('Frequency (THz)')
            ax[2].set_ylabel('Time (ps)')
            ax[4].set_xlabel('LLE Step (sub-sampled)')
            ax[4].set_ylabel('Norm. Comb Pwr')
            ax[6].set_ylabel('Detuning (GHz)')
            ax[0].set_xlim([0,1000])
            [label.set_visible(False) for label in ax[0].get_xticklabels()]
            [label.set_visible(False) for label in ax[2].get_xticklabels()]


            bar_spec = f.colorbar(aa,cax = ax[1],orientation="vertical")
            bar_temp = f.colorbar(bb,cax = ax[3],orientation="vertical")
            bar_spec.set_label('Spec. P (dB)')
            bar_temp.set_label('|E|²')
            tick_locator1 = ticker.MaxNLocator(nbins=4)
            tick_locator2 = ticker.MaxNLocator(nbins=2)
            bar_spec.locator = tick_locator1
            bar_temp.locator = tick_locator2
            bar_spec.update_ticks()
            bar_temp.update_ticks()

            if not pyType == 'jupyter':
                f.show()
            self.fPcomb = f
            self.axPcomb = ax

            return f, ax

        else:
            
            Sspec = go.Heatmap(x=step, y=freq, z=Epb,
                  colorbar=dict(len=0.37, y=0.83, title='Power (dBm)'),
                  yaxis='y3',
                  colorscale='Viridis',
                  zmax = 0,
                  zmin = -120,
                  )

            Stime = go.Heatmap(x = step, y = t, z = E2,
                  colorbar=dict(len=0.34, y=0.47, title = '|E|^2'), 
                   yaxis='y2',
                   colorscale='Viridis')

            Cpwr = go.Scatter(x=step, y=CmbPow,
                              yaxis='y1',
                              name='Comb Power')

            Detuning = go.Scatter(x=step, y=det,
                                  yaxis='y4',
                                  name='Detuning')

            data = [Cpwr, Detuning, Stime, Sspec]
            layout = go.Layout(
                    xaxis=dict(domain=[0, 1],anchor='y1',title= 'LLE Step'),
                    xaxis2=dict(domain=[0, 1],anchor='y1',title= 'couocu', ),
                    xaxis3=dict(domain=[0, 1],anchor='y1'),
                    xaxis4=dict(domain=[0, 1],anchor='y1'),
                    yaxis1=dict(domain=[0, 0.29],title = 'Comb Power',),
                    yaxis2=dict(domain=[0.33, 0.62],anchor='x1',title = 'Fast Time',),
                    yaxis3=dict(domain=[0.66, 1],anchor='x1',title = 'Frequency (THz)',),
                    yaxis4=dict(domain=[0.66, 1],anchor='x1',overlaying='y1',title = 'Detuning (GHz)',side='right'),
                    showlegend=False,)
            fig = go.Figure(data=data, layout=layout)
            iplot(fig)

        self._plotPower = True

    def PlotCombSpectra(self, ind, f=None, ax=None, label=None, pwr='both', do_matplotlib = False):
        '''
        Plot the spectra for a given index in the 1000 sub-sampled LLE steps

        **Input** 

            - ind <ind>: index in the LLE step to plot the spectra
            - f <obj>:  matplotlib figure handle (if None, new figure)
            - ax <obj>: matplotlib axe handle
            - label <str>: label for the legend
            - pwr <str>: 'both', 'ring', 'wg' depending on the spectra wanted (inside the ring, the waveguide or both)

        **Output**

            - freq <numpy.array>: frequency in Hz
            - Sout <numpy.array>: spectral density of power in the waveguide (dBm)
            - Sring <numpy.array>: spectral density of power in the ring (dBm)
            - f <obj>:  matplotlib figure handle
            - ax <obj>: matplotlib axes handle
        '''

        freq = self.sol['freq']*1e-12
        Sring = 30 + 10*np.log10(np.abs(self.sol['Em_probe'][:, ind])**2)
        Sout = 30 + 10*np.log10(np.abs(self.sol['Ewg'][:, ind])**2)
        self.spectra = {'Sout': Sout,
                        'Sres': Sring,
                        'freq': freq*1e-12}


        if not pyType == 'jupyter' or do_matplotlib:
            if f is None and ax is None:
                f, ax = plt.subplots(dpi=120)
            elif f is None and not(ax is None):
                if not type(ax) is list:
                    f = ax.figure
                else:
                    f, ax = plt.subplots(dpi=120)
                    print('Only 1 subplots supported, created a new figure')
            elif not(f is None) and ax is None:
                ax = f.axes[0]
            else:
                if type(ax) is list:
                    f, ax = plt.subplots(dpi=120)

        

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
            if not pyType == 'jupyter':
                f.canvas.draw()
                f.show()
                plt.pause(0.25)

            self.fSpectra = f
            self.axSpectra = ax

            return f, ax
        else:
            trace0 = go.Scatter(x = freq,y = Sring,
                            mode = 'lines',name = 'Res. Power')
            trace1 = go.Scatter(x = freq,y = Sout,
                line = dict(width = 2,dash = 'dash'),
                 mode = 'lines',name = 'Wg. Power')

            if pwr.lower() == 'both':
                data = [trace0, trace1]
            if pwr.lower() == 'ring':
                data = [trace0]
            if pwr.lower() == 'wg':
                data = [trace1]

            layout = dict(xaxis = dict(title = 'Frequency (THz)'),
                  yaxis = dict(title = 'Power (dBm)'),
                  )
            fig = go.Figure(data=data, layout=layout)
            iplot(fig)

        self._plotSpecta = True
        self._indSpectra = ind

    def PlotSolitonTime(self, ind, f=None, ax=None, label=None, do_matplotlib = False):
        '''
        Plot the spectra for a given index in the 1000 sub-sampled LLE step 

        **Input** 

            - ind <ind>: index in the LLE step to plot the spectra
            - f <obj>:  matplotlib figure handle (if None, new figure)
            - ax <obj>: matplotlib axe handle 
            - label <str>: label for the legend

        **Output**

            - τ <obj>: Time in the resonator
            - U <numpy.array>: Temporal Electric field for the given step of the LLE
            - f <obj>: matplotlib figure handle
            - ax <obj>: matplotlib axe handle
        '''

        
        tR = 2*np.pi*self.res['R']*self.res['ng']/self._c0
        freq = self.sol['freq']
        f, ax = plt.subplots(dpi=120)
        τ = np.linspace(-0.5, 0.5, freq.size) * tR
        U = np.abs(self.sol['u_probe'][:,ind])**2

        self.fasttime ={'U': U, 
                        'tau': τ}
        if not pyType == 'jupyter' or do_matplotlib:
            ax.plot(τ*1e12 , U/U.max())
            ax.set_xlabel('Time (ps)')
            ax.set_ylabel('Soliton Energy (a.u)')
            if not pyType == 'jupyter':
                f.show()
            return f, ax
        else:
            trace0 = go.Scatter(x = τ*1e12,y = U,
                            mode = 'lines')
            data = [trace0]
            layout = dict(xaxis = dict(title = 'Time (ps)'),
                  yaxis = dict(title = '|E|^2 (norm)'),
                  )
            fig = go.Figure(data=data, layout=layout)
            iplot(fig)

        self._plotTime = True
        self._indTime = ind

    def SaveResults(self, fname, path='./'):
        '''
        Save the whole class with pickle to be able to easilly call it back or retrieve the results after saving

        **Input**

            - fname <str>: name to save. The '.pkl' extension will be added
            - path <str>: path to save the results (defaults './')
        '''


        to_save = copy(self)
        to_save.sim.pop('domega_disp', None)
        fname = path + fname + '.pkl'
        print(fname)
        pkl.dump(to_save, open(fname,'bw'))

    def SavePlots2Pdf(self,basename):
        if self._plotPower:
            fpwr, axpwr = self.PlotCombPower(do_matplotlib = True)
            Latexify(figname = basename + 'CombPower', fig = fpwr, frmt = 'pdf')
        if self._plotSpecta:
            fspec, axspec = self.PlotCombSpectra(self._indSpectra, do_matplotlib = True)
            Latexify(figname = basename + 'CombSpectra', fig = fspec, frmt = 'pdf')
        if self._plotTime:
            ftime, axtome = self.PlotSolitonTime(self._indTime, do_matplotlib = True)
            Latexify(figname = basename + 'FastTime', fig = ftime, frmt = 'pdf')

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
        
        return to_print
