import numpy as np
from ._analyzedisp import AnalyzeDisp
import scipy.interpolate as itp
import scipy.optimize as optm
import scipy.fftpack as fft
from scipy import constants as cts
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
from prettytable import PrettyTable
from plotly import tools
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
import plotly.graph_objs as go

with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)
    import h5py
import ipdb


backend = mpl.get_backend()
path_juliaScript = os.path.dirname(os.path.abspath(__file__))
path_juliaScript = os.path.join(path_juliaScript, 'ComputeLLE.jl')

print('*'*60)
print('Julia solver:')
print(path_juliaScript)
print('*'*60)



# delete folder as the temp file will not be in it but will be used as a prefix
# print('-'*50)
# print('Path to temporary dir to save .h5 files with prefix:')
# print(tmp_dir)
# print('-'*50)
# print('Path to Julia script: ')
# print(path_juliaScript)
# print('-'*50)


main_tmp_dir = tempfile.mkdtemp()
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

class LLEsolver(object):
    '''
    Class to solve the Lugiato Lefever Equation
    Initialization input ([]=facultative):

    **res <dict>**

        - Qi <float>: intrinsic Q of the resonator
        - Qc <float>: coupling Q of the resonator
        - R <float>: ring radius
        - gamma <float>: Non linear index of the material
        - dispfile <str> : str pointing to a .csv file where the azimuthal mode orders and corresponding resonances are saved

    **sim <dict>**

        - Tscan <float>: length of the simulation (in unit of round trip)
        - mu_fit <list>: number of mode to fit
        - mu_sim <list>: number of mode to simulate
        - domega_init <float>: initial detuning of the pump
        - domega_end <float>: final detuning of the pump
        - [domega_stop] <float>: where to stop the scan in detuning but keep doing the simulation

    **debug <bool>**: Save a trace in a logfile in the working directory of the different actions pyLLE perform (default = False)
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
        self._res = kwargs.get('res', {})
        self._sim = kwargs.get('sim', {})
        # self._sim_norm = kwargs.get('sim_norm', None)
        self._debug = kwargs.get('debug', False)
        self._plotPower = None
        self._plotSpecta = None
        self._indSpectra = 0
        self._plotTime = None
        self._indTime = 0

        # find all the needed parameters
        # ------------------------------------------------------------
        assert 'Qi' in self._res.keys(), 'Please provide Qi'
        assert 'Qc' in self._res.keys(), 'Please provide Qc'
        assert 'R' in self._res.keys(), 'Please provide R'
        assert 'Tscan' in self._sim.keys(), 'Please provide Tscan'
        assert 'dispfile' in self._res.keys(), 'Please provide dispfile'

        if not "D1_manual" in self._sim.keys():
            self._sim["D1_manual"] = None

        # -- check if wavelength or frequency
        # ------------------------------------------------------------
        if not('f_pmp' in self._sim.keys()):
            self._sim['f_pmp'] = self._c0/self._sim['lbd_pmp']
        if not('mu_fit' in self._sim.keys()):
            self._sim['mu_fit'] = [None, None]


        # -- Make sur than each pump has a power --
        if not type(self._sim['f_pmp'])==list:
            self._sim['f_pmp']= [self._sim['f_pmp']]
        if not type(self._sim['Pin'])==list:
            self._sim['Pin']= [self._sim['Pin']]
        assert len(self._sim['Pin']) == len(self._sim['f_pmp']), "The number of *pump* and *pump power* has to be the same"




        # -- Translate the dict from greek sign =--
        # ------------------------------------------------------------
        self._sim = self._Translator(self._sim)
        self._res = self._Translator(self._res)


        # -- Fix missing keys --
        # ------------------------------------------------------------
        if not ('domega_stop' or 'δω_stop') in self._sim.keys():
            try:
                self._sim['domega_stop'] = self._sim['domega_end']

            except:
                self._sim['δω_stop'] = self._sim['δω_end']



        if not np.diff(self._sim['mu_sim'])[0] % 2 == 0:
            print(r"/!\ Simulation modes length need to be odd")
            print('Modification of the simulation modes suchat that:')
            self._sim['mu_sim'] = [self._sim['mu_sim'][0], self._sim['mu_sim'][1]+1]
            print('μ_sim = [{:.0f} {:.0f}]'.format(self._sim['mu_sim'][0],
                                                  self._sim['mu_sim'[1]] ))


        # -- handle the detuning of aux pump--
        # ------------------------------------------------------------
        if not 'domega' in self._sim.keys():
            self._sim['domega']  = [0]*len(self._sim['f_pmp'])
            self._sim['ind_pump_sweep'] = 0
        else:
            domega = self._sim['domega']
            id_sweep = [ii for ii in range(len(domega)) if domega[ii] == None]
            assert len(id_sweep)>=1, 'Need to provide at least None possition in domega to determine which pump is swept'

            self._sim['ind_pump_sweep'] = id_sweep
        assert len(self._sim['f_pmp']) == len(self._sim['domega']), "Please provide the same number of pump fixed detuning than pump frequencies. Detuning = None for the pump which is swept"
        assert len(self._sim['Pin']) == len(self._sim['domega']), "Please provide the same number of pump fixed detuning than pump power. Detuning = None for the pump which is swept"

        if not type(self._sim['ind_pump_sweep'])==list:
            self._sim['ind_pump_sweep'] = [self._sim['ind_pump_sweep']]

        # -- handle the phase of the pump--
        # ------------------------------------------------------------
        if not 'phi_pmp' in self._sim.keys():
            self._sim['phi_pmp']  = [0]*len(self._sim['f_pmp'])
        else:
            if not type(self._sim['phi_pmp'])==list:
                self._sim['phi_pmp'] = [self._sim['phi_pmp']]
        assert len(self._sim['f_pmp']) == len(self._sim['phi_pmp']), "Please provide the same number of pump phase than pump frequencies. Detuning = None for the pump which is swept"



        # -- Setup the Logger ---
        # ------------------------------------------------------------
        if self._debug:
            self._logger = MyLogger("LLE.log")
            self._logger.info('__init__', 'New LLE')
        else:
            self._logger = None

    def _Translator(self,D):
        Dnew = {}
        self._did_translate = False
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
                    self._did_translate = True
                    new_k += self._greek[char]
                else:
                    new_k += char
            Dnew[new_k] = D[k]

        return Dnew

    def Analyze(self, plot=True, display_center = True, label=None, plottype='all', verbose = False):
        '''
        Call pyLLE.analyzedisp.AnalyzeDisp to get the dispersion of the resonator we want to simulate
        '''

        fig = None


        if self._debug:
            Info = '\n\tFilename: {}\n'.format(self._res['dispfile'])
            for ff in self._sim['f_pmp']:
                Info += '\tf_pmp: {:.3f} THz'.format(ff*1e-12)
            self._logger.info('LLEsovler.Analyze', "Analyzing dispersion file" + Info)


        # -- Initialized usufu diserspions data --
        # ------------------------------------------------------------
        self._sim['ind_pmp'] = []
        self._sim['Dint'] = []
        # self._res['ng'] = []
        # self._res['neff'] = []



        # -- Get dispersion for each pump --
        # ------------------------------------------------------------
        # for ff in self._sim['f_pmp']:
        fpmp = copy(self._sim['f_pmp'])
        self._analyze = AnalyzeDisp(file= self._res['dispfile'],
                                     fpmp = fpmp,
                                     rM_fit = self._sim['mu_fit'],
                                     rM_sim = self._sim['mu_sim'],
                                     R = self._res['R'],
                                     D1_manual = self._sim["D1_manual"],
                                     fig = fig,
                                     plottype = plottype,
                                     logger = self._logger,
                                     pyType = pyType)

        self._analyze.GetDint()
        # if verbose:
        #     self._analyze.DisplayParam()
        if plot:
            fig = self._analyze.DisplayPlot(display_center)


        #
        self._sim['Dint'] = self._analyze.Dint_sim[0]
        self._sim['D1'] = self._analyze.D1[0]
        self._res['ng'] = self._analyze.ng_pmp
        self._sim['ind_pmp'] = self._analyze.pmp_ind
        self._sim['f_center'] = self._analyze.fpmp[-1]
        self._sim['f_pmp'] = self._analyze.fpmp[:-1]
        self._sim['mu_sim_center'] = [self._sim['mu_sim'][0]+self._sim['ind_pmp'][0],
                                    self._sim['mu_sim'][1]+self._sim['ind_pmp'][0]]
        self._sim['D1_center'] = self._analyze.D1[-1]
        self._sim['FSR'] = [dd/(2*np.pi) for dd in self._analyze.D1[:-1]]
        self._sim['FSR_center'] = self._analyze.D1[-1]/(2*np.pi)



        if not "DKS_init" in self._sim.keys():
            self._sim['DKS_init']= np.zeros(self._analyze.Dint_sim[0].size)


        disp = {
                'ng': self._analyze.ng,
                'D': self._analyze.D*1e6,
                "D1": self._analyze.D1[0],
                'neff': self._analyze.neff,
                'freq': self._analyze.rf,
                'freq_sim': self._analyze.freq_fit,
                'Dint': self._analyze.Dint[:-1],
                'Dint_sim': self._analyze.Dint_sim[:-1],
                'Dint0': self._analyze.Dint_sim[-1],
                'FSR': self._sim['FSR'],
                'FSR_center': self._sim['FSR_center']
                }

        res = {
            'R': self._res['R'],
            'Qi': self._res['Qi'],
            'Qc': self._res['Qc'],
            'γ': self._res['gamma'],
            'dispfile': self._res['dispfile'],
            }

        sim = {
            "Pin": self._sim['Pin'],
            "Tscan": self._sim['Tscan'],
            "f_pmp": self._sim['f_pmp'],
            "f_center": self._sim['f_center'],
            "δω_init": self._sim['domega_init'],
            "δω_end": self._sim['domega_end'],
            "δω": self._sim['domega'],
            "μ_sim": self._sim['mu_sim'],
            "μ_fit": self._sim['mu_fit'],
            "μcenter": int(-1*self._sim['mu_sim_center'][0]),
            "ind_pmp": self._sim['ind_pmp'],
            "ind_pump_sweep": self._sim['ind_pump_sweep'],
            "phi_pmp": self._sim['phi_pmp'],
            "DKS_init": self._sim['DKS_init'],
        }



        did_translate = not(self._did_translate)

        if not self._did_translate:
            sim = self._Translator(sim)

        self.sim = _dic2struct(sim,which = 'sim', do_tr = did_translate)
        self.disp = _dic2struct(disp,which = 'disp', do_tr = did_translate)
        self.res = _dic2struct(res, which = 'res', do_tr = did_translate)

        return fig

    def Setup(self, tmp_dir = None, verbose = False):
        '''
        Setup the simulation for the Julia back-end.
        Save the two main dictionary self._sim and self._res into a readable hdf5 file for Julia in the temporary location define by the os
        '''
        if tmp_dir == None:
            self.tmp_dir = main_tmp_dir
            # os.rmdir(os.path.join(tmp_dir')
        else:
            self.tmp_dir = tmp_dir

        def CleanUpHDF5():

            try:
                os.remove(self.tmp_dir + 'ParamLLEJulia.h5')
            except:
                pass

            try:
                os.remove(self.tmp_dir + 'ResultsJulia.h5')
            except:
                pass
            try:
                os.remove(self.tmp_dir + 'log.log')
            except:
                pass

        def LoggerDic():
            dic_sim = {'Pin': ('Pin',1e3, 'mW'),
                        'Tscan': ('Tscan',1e-6, 'x1e6 Round Trip'),
                        'domega_init': (u'\u03B4\u03C9_init',1e-9/(2*np.pi), u'x2\u03C0 GHz'),
                        'domega_end': (u'\u03B4\u03C9_end',1e-9/(2*np.pi), u'x2\u03C0 GHz'),
                        'domega_stop': (u'\u03B4\u03C9_stop',1e-9/(2*np.pi), u'x2\u03C0 GHz'),
                        'f_pmp': ('f_pmp',1e-12, 'THz'),
                        'domega_aux': (u'\u03B4\u03C9_aux',1e-9/(2*np.pi), u'x2\u03C0 GHz'),
                        'ind_pump_sweep': ('ind_pump_sweep',1, ''),
                        'mu_sim': (u'\u03BC_sim',1, ''),
                        'mu_fit': (u'\u03BC_fit',1, ''),}

            try:
                if len(self._res['Qc']) >1:
                    dic_res = {'R': ('R',1e6, 'µm'),
                                'Qi': ('Qi',1e-6, 'M'),
                                'gamma': (u'\u03B3', 1, ''),}
                else:
                    dic_res = {'R': ('R',1e6, 'µm'),
                                'Qi': ('Qi',1e-6, 'M'),
                                'Qc': ('Qc',1e-6, 'M'),
                                'gamma': (u'\u03B3', 1, ''),}
            except:
                dic_res = {'R': ('R',1e6, 'µm'),
                            'Qi': ('Qi',1e-6, 'M'),
                            'Qc': ('Qc',1e-6, 'M'),
                            'gamma': (u'\u03B3', 1, ''),}

            return dic_sim, dic_res

        def PrintLogger(dic_sim, dic_res):
            self._sim['debug'] = int(self._debug)
            Info = '-- Solving standard LLE --\n'
            Info += '\tSimulation Parameters\n'
            for k, it in self._res.items():
                if k in dic_res.keys():
                    if it == None:
                        Info +='\t\t{} = {}\n'.format(dic_res[k][0], 'None')
                    else:
                        Info +='\t\t{} = {:.2f} {}\n'.format(dic_res[k][0], it*dic_res[k][1],dic_res[k][2])


            Info += '\tSimulation Parameters\n'
            for k, it in self._sim.items():
                if k in dic_sim.keys():
                    if type(it) is list:
                        if dic_sim[k][0] == u'\u03BC_sim' or dic_sim[k][0] == u'\u03BC_fit':
                            if it == [None, None]:
                                Info += '\t\t{} = [None,None] {}\n'.format(dic_sim[k][0],
                                                                    dic_sim[k][2])
                            else:
                                Info += '\t\t{} = [{:.2f},{:.2f}] {}\n'.format(dic_sim[k][0],
                                                                    it[0]*dic_sim[k][1],
                                                                    it[1]*dic_sim[k][1],
                                                                    dic_sim[k][2])
                        else:
                            cnt = 0
                            for iitt in it:
                                if iitt == None:
                                    Info += '\t\t{}[{}] = {}\n'.format(dic_sim[k][0],
                                                                        cnt,
                                                                       'None')
                                else:
                                    Info += '\t\t{}[{}] = {:.2f} {}\n'.format(dic_sim[k][0],
                                                                            cnt,
                                                                            iitt*dic_sim[k][1],
                                                                            dic_sim[k][2])
                                cnt += 1
                    else:
                        if it == None:
                            Info +='\t\t{} = {} {}\n'.format(dic_sim[k][0], 'None',dic_sim[k][2])
                        else:
                            Info +='\t\t{} = {:.2f} {}\n'.format(dic_sim[k][0], it*dic_sim[k][1],dic_sim[k][2])

            print(Info)
            if self._debug:
                try:
                    self._logger.info('LLEsovler.Setup', Info)
                except:
                    Info = ''.join([self._greek[ii] if ii in self._greek.keys() else ii for ii in Info])
                    self._logger.info('LLEsovler.Setup', Info)

        def SetupHDF5():
            # -- create h5file --
            if verbose:
                print('HDF5 parameter file can be foud in: {}'.format(self.tmp_dir + 'ParamLLEJulia.h5'))
            h5f = h5py.File(self.tmp_dir + 'ParamLLEJulia.h5', 'w')
            if self._debug:
                self._logger.info('LLEsovler.Setup','Saving parameters in: {}'.format(self.tmp_dir + 'ParamLLEJulia.h5'))


            h5f.create_group('sim')
            h5f.create_group('res')

            for key, it in self._sim.items():
                if not key == 'domega_disp':
                    if type(it) is str:
                        it = np.string_(it)
                    if type(it) is list:
                        try:
                            if None in it:
                                it = [0 if iitt is None else iitt for iitt in it]
                        except:
                            pass
                    else:

                        try:
                            if it == None:
                                it = 0
                        except:
                            pass
                    dset = key
                    h5f.create_dataset(f'sim/{dset}', data=[it])
                if key == "DKS_init":

                    dset = 'DKSinit_real'
                    h5f.create_dataset(f'sim/{dset}', data=np.real(it))
                    dset = 'DKSinit_imag'
                    h5f.create_dataset(f'sim/{dset}', data=np.imag(it))




            for key, it in self._res.items():
                if not key == 'domega_disp':
                    if type(it) is str:
                        it = np.string_(it)
                    if type(it) is list:
                        if None in it:
                            it = [0 if iitt is None else iitt for iitt in it]
                    else:
                        if it == None:
                            it = 0
                    h5f.create_dataset('res/{}'.format(key), data=[it])
            h5f.close()



        # delay a bit the start
        time.sleep(0.5)
        CleanUpHDF5()
        # make sure the cleanup happened
        time.sleep(1)
        dic_sim, dic_res = LoggerDic()
        if verbose or self._debug:
            PrintLogger(dic_sim, dic_res)
        SetupHDF5()

    def SolveTemporal(self, verbose = False, tol = 1e-3, maxiter = 6, step_factor = 0.1):
        '''
        Call Julia to solve the LLE
        '''
        self._solver = 'temporal'


        def DisplaySim():

            if self._debug:
                self._logger.info('LLEsovler.SolveTemporal','Solving Temporal LLE with Julia....')
                Info = 'tol = {} -- maxiter = {} step_factpr = {}'.format(tol, maxiter, step_factor)
                self._logger.info('LLEsovler.SolveTemporal',Info)


            hLine = '-'*70

            print(hLine)

            date = time.localtime()
            date = "{}-{:0>2}-{:0>2} ".format(date.tm_year,
                                              date.tm_mon,
                                              date.tm_mday) +\
                "{:0>2}:{:0>2}:{:0>2}".format(date.tm_hour,
                                              date.tm_min,
                                              date.tm_sec)

            print(date)

        def LaunchJulia():

            julia = 'julia'

            try:
                self.JuliaSolver = sub.Popen(julia, stdout=sub.PIPE, stderr=sub.PIPE)
            except:
                raise ValueError('julia is not installed on the system path. ',
                                 'Please add julia to the path or re-install ',
                                 'and check the option to add to path.')

            command = [julia, path_juliaScript , self.tmp_dir, str(tol), str(maxiter), str(step_factor)]
            self.JuliaSolver = sub.Popen(command, stdout=sub.PIPE, stderr=sub.PIPE)
            fname = self.tmp_dir + 'log.log'
            if verbose:
                print('Launching Julia....')
                print('Temp file can be found in: {}'.format(fname))

            start = time.time()
            return fname, start

        def Pbar(perc, pgrs, tb_up, conv_err):
            bar_width = 50
            pgrs = '*'*int(np.floor(perc/2))
            width = ' '* (bar_width - len(pgrs))
            perc_str = ' {}%'.format(perc)
            line = 'Computing LLE [' + pgrs  + width  + ']' +  perc_str
            if conv_err:
                line = line + ' /!\ Convergence issue'
                if self._debug:
                    self._logger.info('LLEsovler.SolveTemporal','/!\ Convergence issue')
            length = len(line)
            return line, length, pgrs, tb_up

        def ProbeProgress(fname, start):
            tb_up = 2
            pgrs = ''
            perc_old = -1
            perc = -1
            conv_err = False
            line = ''
            timenow = time.time()
            # wait for the solver to actually start
            while not perc == 0 and self.JuliaSolver.poll() == None:
                try:
                    perc = int(open(fname).readlines()[-1].strip())
                except Exception as e:
                    pass
                    # print('line 654')
                    # print(e)

            #fetch if any errors:
            err = False
            if not self.JuliaSolver.poll() == None:
                _, err = self.JuliaSolver.communicate()
                if not err.decode() == '':
                    print('!!! JULIA ERROR !!!')
                    print(err.decode())
                    err = True
                print("Error: {}".format(err))

            if not err:
                if verbose:
                    print('Launching Julia: Done')
                perc = -1
                while not perc == 100 and self.JuliaSolver.poll() == None:
                    try:
                        ll = open(fname).readlines()[-1].strip()
                        try:
                            perc = int(ll)
                            if not perc_old == perc:
                                line, length, pgrs, tb_up= Pbar(perc, pgrs, tb_up, conv_err)
                                print('\r' + line, end = '')
                                perc_old = perc

                        except Exception as e:
                            print('line 681')
                            print(e)
                            if ll.split()[0] == 'Failed':
                                conv_err = True

                    except Exception as e:
                        print('line 686')
                        print(e)
                        pass
                    time.sleep(1)
                line, length, pgrs, tb_up= Pbar(100, pgrs, tb_up, conv_err)
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

        DisplaySim()
        time.sleep(2)
        JuliaLog, start_time = LaunchJulia()
        ProbeProgress(JuliaLog, start_time)

    def SolveSteadyState(self, do_plot = True):
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
        μ_sim = copy(self._sim['mu_sim'])
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
        Pin = self._sim['Pin']
        γ = self._res['gamma']
        L = 2*np.pi*self._res['R']
        ω0 = self._sim['f_pmp'][0]*2*np.pi
        Q0 = self._res['Qi']
        Qc = self._res['Qc']
        tR = L*self._res['ng']/self._c0
        α = 1/2 * (ω0/Q0 + ω0/Qc) * tR
        θ = ω0/Qc*tR
        δω = -self._sim['domega']* tR

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
        ut0 = np.sqrt(α/γ/L) * (Φ0+ B0 * np.exp(1j*φ0) * sech(0.5*B0*np.sqrt(α/(np.abs(β2)*L))*τ))
        Em0 = fft.fftshift(fft.fft(ut0))
        x_init = np.concatenate((Em0.real, -Em0.imag))

        # -- Define the Steady State LLE Equation --
        φ = -α + 1j*δω - 1j*self._sim["dphi"]
        Em= lambda xx:  xx[0:int(xx.shape[0]/2)] + 1j*xx[int(xx.shape[0]/2)]
        Ut=  lambda xx:  fft.ifft(Em(xx));
        fm= lambda xx: φ*Em(xx) + nlc*fft.fft(np.abs(Ut(xx))**2*Ut(xx)) + Ein_couple;
        fvec= lambda xx: np.concatenate((fm(xx).real, fm(xx).imag))

        # -- Solver the Steady State --
        out = optm.root(fvec, x_init, method='lm', jac=None, tol=1e-8)
        Ering = Em(out.x)/μ.size
        Ewg = Ein/μ.size -Ering*np.sqrt(θ)

        self.steady = {'Ering':Ering, 'Ewg':Ewg, 'Uring': fft.ifft(Ering), 'freq': ν, 'tau': τ}
        if not pyType == 'jupyter':
            if do_plot:
                f, ax = plt.subplots(dpi=120)
                ax.plot(1e-12*ν, 30 + 10*np.log10(np.abs(Ewg)**2), label='Waveguide')
                ax.plot(1e-12*ν, 30 + 10*np.log10(np.abs(Ering)**2), label='Ring')
                ax.legend()
                f.show()
                return f, ax

        else:
            if do_plot:
                trace0 = go.Scatter(x = 1e-12*ν,y = 30 + 10*np.log10(np.abs(Ering)**2),
                                mode = 'lines', name='Res. Power')
                trace1 = go.Scatter(x = 1e-12*ν,y = 30 + 10*np.log10(np.abs(Ewg)**2),
                                mode = 'lines', name='Out Power')
                data = [trace1, trace0]
                layout = dict(xaxis = dict(title = 'Frequency (THz)'),
                      yaxis = dict(title = 'Power (dBm)'),
                      )
                fig = go.Figure(data=data, layout=layout)
                iplot(fig)
                return fig

    def RetrieveData(self):
        '''
        Load the output hdf5 saved by julia and transform it in a user-friendly dictionary to be more pythonistic
        '''

        time.sleep(2)
        drct = self.tmp_dir
        with h5py.File(self.tmp_dir + 'ResultsJulia.h5', 'r') as S:
            if self._debug:
                self._logger.info('LLEsovler.RetrieveData','Retrieving results from Julia in {}'.format(self.tmp_dir + 'ResultsJulia.h5'))
            self._sol = {}
            keys = ['u_probe', 'driving_force']
            for k in keys:
                rl = 'Results/{}Real'.format(k)
                im = 'Results/{}Imag'.format(k)
                self._sol[k] = S[rl][:] + 1j*S[im][:]

            self._sol['detuning'] = S["Results/detuningReal"][:]
            self._sol['time'] = S["Results/t_simReal"][:]
            self._sol['kappa_ext'] = S["Results/kappa_extReal"][:]
            self._sol['kappa0'] = S["Results/kappa_0Real"][:]


        os.remove(self.tmp_dir + 'ParamLLEJulia.h5')
        os.remove(self.tmp_dir + 'ResultsJulia.h5')



        # Now Compute the uselfull parameters
        # ---------------------------------------------
        ind_pump_sweep = self._sim['ind_pump_sweep']
        μcent = -self._sim['mu_sim_center'][0]
        # ind_pmp = [self._sim['ind_pmp'][ii] for ii in ind_pump_sweep]
        det = 1e-9*self._sol['detuning']/(2*np.pi)


        L = self._res['R']
        ng = self._res['ng']
        Qc = self._res['Qc']
        Q0 = self._res['Qi']
        tR = L*ng/cts.c
        FSR = cts.c/(ng*L)
        μlength =  self._sim['mu_sim_center'][1]*2+1

        κext = self._sol['kappa_ext']
        κ0 = self._sol['kappa0']
        α = κext+κ0

        u_probe = self._sol['u_probe']
        driving_force = self._sol['driving_force']
        D1 = self._sim['D1_center']

        ωcenter = self._sim['f_center']*2*np.pi
        μcenter = np.linspace(self._sim['mu_sim_center'][0],
                              self._sim['mu_sim_center'][1],
                              μlength)

        u_probe = self._sol['u_probe']
        Awg = 1j*np.zeros(u_probe.shape)
        Ewg = 1j*np.zeros(u_probe.shape)
        Ecav = 1j*np.zeros(u_probe.shape)
        Acav = 1j*np.zeros(u_probe.shape)
        Pcomb = np.zeros(u_probe.shape[1])
        for ii in range(u_probe.shape[1]):
            wg = driving_force[:,ii]*np.sqrt(1-κext)
            cav = np.sqrt(κext)*u_probe[:,ii]*np.exp(1j*np.pi)
            Awg[:, ii] = (wg + cav)/np.sqrt(μlength)
            Acav[:, ii] = np.sqrt(α[0]/2)*u_probe[:,ii]*np.exp(1j*np.pi)/np.sqrt(μlength)
            Ewg[:, ii] = np.fft.fftshift(np.fft.fft(Awg[:, ii]))/np.sqrt(μlength)+1e-12
            Ecav[:, ii] = np.fft.fftshift(np.fft.fft(Acav[:, ii]))/np.sqrt(μlength)
            Pcomb[ii] = np.sum(np.abs(Ecav[:,ii])**2,0) - np.sum([np.abs(Ecav[μcent+jj,ii])**2 for  jj in np.unique(self._sim['ind_pmp'])])
        Pwg = np.sum(np.abs(Awg)**2,0)
        Pcav = np.sum(np.abs(Acav)**2,0)

        # -- Comb Frequency --
        μ = self._sim['mu_sim']
        μ = np.arange(μ[0], μ[1]+1)
        FSR = self._sim['FSR']
        self._sol['freq'] = self._sim['f_pmp'][0]+μ*FSR[0]

        sol = {'δfreq': self._sol['detuning']/(2*np.pi),
               'time': self._sol['time'],
               'Awg': Awg,
               'Acav': Acav,
               'Ewg': Ewg,
               'Ecav': Ecav,
               'Pcomb': Pcomb,
               'Pwg': Pwg,
               'Pcav': Pcav,
               'μ':μ,
               'freq': self._sol['freq']
               }

        did_translate = not(self._did_translate)
        if not self._did_translate:
            sol = self._Translator(sol)
        self.sol = _dic2struct(sol, which = 'sol', do_tr = did_translate)

    def PlotSpectraMap(self, do_matplotlib = False):

        tr = [go.Heatmap(z = 10*np.log10(np.abs(self.sol.Ewg)**2)+30,
                 zmax = 0,
                 zmin = -90,
                colorbar = dict(title = 'In-wg Power (dBm)'))]
        # tr = [ go.Scatter(y = np.abs(wg)**2)]
        fig = go.FigureWidget(data = tr)
        _ = fig.update_xaxes(title = 'LLE step')
        _ = fig.update_yaxes(title = 'Frequency (THz)')
        # fig.show()

        return fig

    def PlotTimeMap(self, do_matplotlib = False):
        pass

    def PlotCombPower(self, do_matplotlib = False, which = 'all', xaxis = 'steps' ):
        '''

        '''
        tr = []
        if xaxis.lower() == 'steps':
            x = np.linspace(0,999,1000)
            xlabel = 'LLE steps'
        elif xaxis.lower() == 'detuning':
            x = 1e-9*self._sol['detuning']/(2*np.pi)
            xlabel = 'Detuning (GHz)'

        if not pyType == 'jupyter' or do_matplotlib:
            fig, ax = plt.subplots()
            if which.lower() == 'all':
                ax.plot(x, self.sol.Pcav)
                ax.plot(x, self.sol.Pwg)
                ax.plot(x, self.sol.Pcomb)

            if which.lower() == 'comb':
                ax.plot(x, self.sol.Pcomb)
            if which.lower() == 'waveguide':
                ax.plot(x, self.sol.Pwg)
            if which.lower() == 'cavity':
                ax.plot(x, self.sol.Pcav)

            ax.legend()
            ax.set_xlabel(xlabel)
            ax.set_ylabel('Power (dBm)')

        else:
            if which.lower() == 'all':
                tr += [go.Scatter(x = x, y = self.sol.Pcav, name = 'Pcav')]
                tr += [go.Scatter(x = x, y = self.sol.Pwg, name = 'Pwg')]
                tr += [go.Scatter(x = x, y = self.sol.Pcomb, name = 'Pcomb')]

            if which.lower() == 'comb':
                tr += [go.Scatter(x = x, y = self.sol.Pcomb, name = 'Pcomb')]
            if which.lower() == 'waveguide':
                tr += [go.Scatter(x = x, y = self.sol.Pwg, name = 'Pwg')]
            if which.lower() == 'cavity':
                tr += [go.Scatter(x = x, y = self.sol.Pcav, name = 'Pcav')]


            fig = go.FigureWidget(data = tr)
            _ = fig.update_xaxes(title = xlabel)
            _ = fig.update_yaxes(title = 'Power (W)')


        return fig

    def PlotCombSpectra(self, ind, do_matplotlib = False, plot = True, style = 'comb', xaxis = 'freq', where = 'waveguide', floor = -100):
        # freq = self.sol['freq']*1e-12
        # Sring = 30 + 10*np.log10(np.abs(self.sol['Em_probe'][:, ind])**2)
        # Sout = 30 + 10*np.log10(np.abs(self.sol['Ewg'][:, ind])**2)
        # self.spectra = {'Sout': Sout,
        #                 'Sres': Sring,
        #                 'freq': freq*1e-12}
        if xaxis.lower() == 'frequencies' or xaxis.lower() == 'freq':
            x = self.sol.freq*1e-12
            xlabel = 'Frequency (THz)'
        elif xaxis.lower() == 'modes':
            x = np.arange(self._sim['mu_sim'][0],self._sim['mu_sim'][-1])  - self._sim['ind_pmp'][0]
            xlabel = 'Mode Numbers'

        if where == 'waveguide' or where =='wg':
            y = 10*np.log10(np.abs(self.sol.Ewg[:,ind])**2)+30
            name = ['Waveguide']
        if where == 'cavity' or where =='cav':
            y = 10*np.log10(np.abs(self.sol.Ecav[:,ind])**2)+30
            name = ['Cavity']
        if where == 'both':
            y = 10*np.log10(np.abs(self.sol.Ecav[:,ind])**2)+30
            y2 = 10*np.log10(np.abs(self.sol.Ewg[:,ind])**2)+30
            name = ['Cavity', 'Waveguide']

        if style == 'comb':
            x_ = np.zeros(x.size*3)
            x_[::3] = x
            x_[1::3] = x
            x_[2::3] = x
            x = x_

            y_ = np.zeros(y.size*3)
            y_[::3] = floor
            y_[1::3] = y
            y_[2::3] = floor
            y = y_

            if where == 'both':
                y2_ = np.zeros(y2.size*3)
                y2_[::3] = floor
                y2_[1::3] = y2
                y2_[2::3] = floor
                y2 = y2_


        if not pyType == 'jupyter' or do_matplotlib:
            fig, ax = plt.subplots()
            ax.plot(x, y, label = name[0])
            if where == 'both':
                ax.plot(x, y2, label = name[1])
            ax.set_xlabel(xlabel)
            ax.set_ylabel('Power (dBm)')
        else:

            tr = [go.Scatter(x=x, y = y,
                            name = name[0])]
            if where == 'both':
                tr += [go.Scatter(x=x, y = y2,
                                name = name[1])]

            fig = go.FigureWidget(data = tr)

            _ = fig.update_xaxes(title = xlabel)
            _ = fig.update_yaxes(title = 'Power (dBm)')

        self._plotSpecta = True
        self._indSpectra = ind



        return fig

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


        tR = 2*np.pi*self._res['R']*self._res['ng']/self._c0
        freq = self.sol['freq']

        τ = np.linspace(-0.5, 0.5, freq.size) * tR
        U = np.abs(self.sol['Acav'][:,ind])**2

        self.fasttime ={'U': U,
                        'tau': τ}
        if not pyType == 'jupyter' or do_matplotlib:
            f, ax = plt.subplots(dpi=120)
            ax.plot(τ*1e12 , U/U.max())
            ax.set_xlabel('Time (ps)')
            ax.set_ylabel('Soliton Energy (a.u)')
            if not pyType == 'jupyter':
                f.show()
            fig = f
        else:
            trace0 = go.Scatter(x = τ*1e12,y = U,
                            mode = 'lines')
            data = [trace0]
            layout = dict(xaxis = dict(title = 'Time (ps)'),
                  yaxis = dict(title = '|E|^2 (norm)'),
                  )
            fig = go.FigureWidget(data=data, layout=layout)


        self._plotTime = True
        self._indTime = ind
        return fig

    def SaveResults(self, fname, path='./'):
        '''
        Save the whole class with pickle to be able to easilly call it back or retrieve the results after saving

        **Input**

            - fname <str>: name to save. The '.pkl' extension will be added
            - path <str>: path to save the results (defaults './')
        '''
        to_save = copy(self)
        # to_save._sim.pop('domega_disp', None)
        # to_save.sim.pop('domega_disp', None)
        del to_save.JuliaSolver
        fname = path + fname + '.pkl'
        print(fname)
        pkl.dump(to_save, open(fname,'bw'))

    def SavePlots2File(self,basename = './', format = 'pdf'):
        if self._plotPower:
            fpwr, axpwr = self.PlotCombPower(do_matplotlib = True)
            Latexify(figname = basename + 'CombPower', fig = fpwr, frmt = format)
        if self._plotSpecta:
            fspec, axspec = self.PlotCombSpectra(self._indSpectra, do_matplotlib = True)
            Latexify(figname = basename + 'CombSpectra', fig = fspec, frmt = format)
        if self._plotTime:
            ftime, axtome = self.PlotSolitonTime(self._indTime, do_matplotlib = True)
            Latexify(figname = basename + 'FastTime', fig = ftime, frmt = format)

    def __repr__(self):
        to_print = ''
        to_print = 'Dispersion from: {}\n\n'.format(self._res['dispfile'])
        to_print += 'Resonator Parameters:\n'
        res_table = PrettyTable(['Parameters', 'Value', 'Units'])
        res_table.add_row(['R', "{:.3f}".format(self._res['R']*1e6),'µm'])
        res_table.add_row(['Qi', "{:.3f}".format(self._res['Qi']*1e-6),'x1e6'])
        res_table.add_row(['Qc', "{:.3f}".format(self._res['Qc']*1e-6),'x1e6'])
        if 'gamma' in self._res:
            res_table.add_row(['γ', "{:.3f}".format(self._res['gamma']),''])
        if 'n2' in self._res:
            res_table.add_row(['n2', "{:.3f}".format(self._res['n2']*1e19),'x1e-19 m2/W'])
        to_print += res_table.get_string()
        to_print += '\n'
        to_print += '\n'

        to_print += 'Simulation Parameters:\n'
        sim_table = PrettyTable(['Parameters', 'Value', 'Units'])
        if 'Pin' in self._sim:
            for pp in self._sim['Pin']:
                sim_table.add_row(['Pin',"{:.3f}".format(pp*1e3),'mW'])
        if 'f_pmp' in self._sim:
            for ff in self._sim['f_pmp']:
                sim_table.add_row(['f_pmp',"{:.3f}".format(ff*1e-12),'THz'])
        if 'μ_sim' in self._sim:
            sim_table.add_row(['μ_sim',"{}".format(self._sim['mu_sim']),''])
        if 'Tscan' in self._sim:
            sim_table.add_row(['Tscan',"{:.3f}".format(self._sim['Tscan']*1e-5),'x1e5 tR'])
        if 'domega_init' in self._sim:
            sim_table.add_row(['δω_init',"{:.3f}".format(self._sim['domega_init']*1e-9),'GHz'])
        if 'domega_end' in self._sim:
            sim_table.add_row(['δω_end',"{:.3f}".format(self._sim['domega_end']*1e-9),'GHz'])
        to_print += sim_table.get_string()
        to_print += '\n'

        return to_print




class _dic2struct():
    def __init__(self, d, which='sim', do_tr=True):
        self._dic = d
        self._which = which
        self._do_tr= do_tr
        self._test  = LLEsolver
        for a, b in d.items():
           setattr(self, a, _dic2struct(b) if isinstance(b, dict) else b)

    def reprRes(self):
        greek ={'α': 'alpha', 'β':'beta','γ': 'gamma',
                    'ω': 'omega','δ': 'd', 'Δ': 'D', 'μ': 'mu',
                    'λ': 'lambda','ε': 'epsilon','φ': 'phi'}
        table = PrettyTable(['Parameter', 'Description', 'Values', 'Units'])
        table.add_row(['R', self._dic['R']*1e6, 'Ring radius', 'µm'])
        Qi = self._dic['Qi']
        exp = np.floor(np.log10(Qi))
        arg = Qi*10**(-1*exp)
        Qi = '{:.3f} x10^{:.0f}'.format(arg, exp)
        table.add_row(['Qi', Qi, 'Intrinsic quality factor', ''])

        Qc = self._dic['Qc']
        exp = np.floor(np.log10(Qc))
        arg = Qc*10**(-1*exp)
        Qc = '{:.3f} x10^{:.0f}'.format(arg, exp)
        table.add_row(['Qc', Qc, 'Coupling quality factor', ''])
        try:
            table.add_row(['γ',self._dic['gamma'],  'Non-linear coefficient', 'W^-1 m^-1'])
        except:
            table.add_row(['γ',self._dic['γ'],  'Non-linear coefficient', 'W^-1 m^-1'])
        table.add_row(['dispfile', self._dic['dispfile'], 'mode/resonance frequency file', ''])
        table.vrules = False
        table.header = True
        table.align = "l"
        table.padding_width = 2
        table.padding_height = 25

        str_table = table.get_string()
        if self._do_tr:
            for ii in greek.items():
                str_table = str_table.replace(ii[0], ii[1] )
        return str_table

    def reprSim(self):
        greek ={'α': 'alpha', 'β':'beta','γ': 'gamma',
                'ω': 'omega','δ': 'd', 'Δ': 'D', 'μ': 'mu',
                'λ': 'lambda','ε': 'epsilon','φ': 'phi'}

        table = PrettyTable(['Parameter', 'Description', 'Values', 'Units'])
        Pin = ['{:.3f}'.format(ii*1e3) for ii in self._dic['Pin']]
        Pin = '[' + ', '.join(Pin) + ']'
        table.add_row(["Pin",Pin + ' (mW)', "Pump power", "W"])
        table.add_row(["Tscan",self._dic['Tscan'], "Simulation time length", "unit of round trip"])

        f_pmp = ['{:.3f}'.format(ii*1e-12) for ii in self._dic['f_pmp']]
        f_pmp = '[' + ', '.join(f_pmp) + ']'
        table.add_row(["f_pmp",f_pmp + ' (THz)',  "Pump frequencies", "Hz"])
        table.add_row(["φ_pmp",self._dic['phi_pmp'], "Pump phase ", "rad"])

        f_center = '{:.3f}'.format(self._dic['f_center']*1e-12)
        table.add_row(["f_center",f_center + ' (THz)', "Center of the sim domain", "Hz"])

        try:
            δω_init = '{:.3f}'.format(self._dic['δω_init']*1e-9/(2*np.pi))
            δω_end = '{:.3f}'.format(self._dic['δω_end']*1e-9/(2*np.pi))
        except:
            δω_init = '{:.3f}'.format(self._dic['domega_init']*1e-9/(2*np.pi))
            δω_end = '{:.3f}'.format(self._dic['domega_end']*1e-9/(2*np.pi))
        δω =[]
        try:
            for ii in self._dic['δω']:
                if ii:
                    δω += ['{:.3f}'.format(ii*1e-9/(2*np.pi))]
                else:
                    δω += ['None']
        except:
            for ii in self._dic['domega']:
                if ii:
                    δω += ['{:.3f}'.format(ii*1e-9/(2*np.pi))]
                else:
                    δω += ['None']
        δω = '[' + ', '.join(δω) + ']'
        table.add_row(["δω_init",δω_init + ' (GHz)', "Start of the frequency sweep", "x2π Hz"])
        table.add_row(["δω_end",δω_end + ' (GHz)', "End of the frequency sweep", "x2π Hz"])
        table.add_row(["δω",δω + ' (GHz)', "Fixed detuning", "x2π Hz"])


        try:
            table.add_row(["μ_sim",self._dic['μ_sim'], "Simulation mode domain", ""])
            table.add_row(["μ_fit",self._dic['μ_fit'], "Fit mode domain", ""])
            table.add_row(["μcenter",self._dic['μcenter'], "Index of the center of the sim domain", ""])
        except:
            table.add_row(["μ_sim",self._dic['mu_sim'], "Simulation mode domain", ""])
            table.add_row(["μ_fit",self._dic['mu_fit'], "Fit mode domain", ""])
            table.add_row(["μcenter",self._dic['mucenter'], "Index of the center of the sim domain", ""])
        table.add_row(["ind_pmp",self._dic['ind_pmp'], "Pump index relative to center of domain", ""])
        table.add_row(["ind_pump_sweep",self._dic['ind_pump_sweep'][0], "Pump index to sweep", ""])
        table.vrules = False
        table.header = True
        table.align = "l"
        table.padding_width = 2
        table.padding_height = 25
        str_table = table.get_string()
        if self._do_tr:
            for ii in greek.items():
                str_table = str_table.replace(ii[0], ii[1] )
        return str_table

    def reprSol(self):
        greek ={'α': 'alpha', 'β':'beta','γ': 'gamma',
                    'ω': 'omega','δ': 'd', 'Δ': 'D', 'μ': 'mu',
                    'λ': 'lambda','ε': 'epsilon','φ': 'phi'}
        table = PrettyTable(['Parameter', 'Description', 'Values', 'Units'])

        try:
            δfreq  = '[{:.3f} ... {:.3f}]'.format(self._dic['δfreq'][0]*1e-9, self._dic['δfreq'][-1]*1e-9)
        except:
            δfreq  = '[{:.3f} ... {:.3f}]'.format(self._dic['dfreq'][0]*1e-9, self._dic['dfreq'][-1]*1e-9)

        table.add_row(["δfreq", δfreq + ' (GHz)',  "Pump detuning", "Hz"])

        time = '[{:.3f} ... {:.3f}]'.format(self._dic['time'][0]*1e6, self._dic['time'][-1]*1e6)
        table.add_row(["time", time + ' (μs)', "Simualtion time", "s"])

        N = self._dic['Awg'].shape
        N = '[{} fast time x {} slow time]'.format(N[0], N[1])
        table.add_row(["Awg", N ,"E. field in time domain in wg,", "V/m"])
        table.add_row(["Acav", N, "E. field in time domain in cav.", "V/m"])
        N = self._dic['Ewg'].shape
        N = '[{} spectra x {} slow time]'.format(N[0], N[1])
        table.add_row(["Ewg",N, "E. field in freq. domain in wg.", "V/m"])
        table.add_row(["Ecav",N, "E. field in freq. domain in cav.", "V/m"])

        Pcomb = '[{:.3f} ... {:.3f}]'.format(self._dic['Pcomb'][0]*1e3, self._dic['Pcomb'][-1]*1e3)
        table.add_row(["Pcomb", Pcomb + ' (mW)', "Intra-cavity comb power", "W"])

        Pwg = '[{:.3f} ... {:.3f}]'.format(self._dic['Pwg'][0]*1e3, self._dic['Pwg'][-1]*1e3)
        table.add_row(["Pwg", Pwg + ' (mW)', "In-waveguide power", "W"])

        Pcav = '[{:.3f} ... {:.3f}]'.format(self._dic['Pcav'][0]*1e3, self._dic['Pcav'][-1]*1e3)
        table.add_row(["Pcav", Pcav + ' (mW)', "Intracavity power", "W"])

        try:
            μ = '[{:.0f} ... {:.0f}]'.format(self._dic['μ'][0], self._dic['μ'][-1])
        except:
            μ = '[{:.0f} ... {:.0f}]'.format(self._dic['mu'][0], self._dic['mu'][-1])
        table.add_row(["μ", μ, "Frequency comb modes index", ""])

        freq = '[{:.3f} ... {:.3f}]'.format(self._dic['freq'][0]*1e-12, self._dic['freq'][-1]*1e-12)
        table.add_row(["freq", freq + ' (THz)', "Frequency comb frequencies", "Hz"])
        table.vrules = False
        table.header = True
        table.align = "l"
        table.padding_width = 2
        table.padding_height = 25
        str_table = table.get_string()
        if self._do_tr:
            for ii in greek.items():
                str_table = str_table.replace(ii[0], ii[1] )
        return str_table

    def reprDisp(self):
        greek ={'α': 'alpha', 'β':'beta','γ': 'gamma',
                    'ω': 'omega','δ': 'd', 'Δ': 'D', 'μ': 'mu',
                    'λ': 'lambda','ε': 'epsilon','φ': 'phi'}
        table = PrettyTable(['Parameter', 'Description', 'Values', 'Units'])
        ng = '[{:.3f} ... {:.3f}]'.format(self._dic['ng'][0], self._dic['ng'][-1])
        table.add_row(['ng', ng, 'Group Index', ''])
        D  = '[{:.3f} ... {:.3f}]'.format(self._dic['D'][2], self._dic['D'][-3])
        table.add_row(['D', D, 'Dispersion', '(ps/nm/km)'])

        neff  = '[{:.3f} ... {:.3f}]'.format(self._dic['neff'][1], self._dic['neff'][-2])
        table.add_row(['neff', neff, 'Effctive Index', ''])
        freq  = '[{:.3f} ... {:.3f}]'.format(self._dic['freq'][1]*1e-12,  self._dic['freq'][-2]*1e-12)
        table.add_row(['freq', freq + ' (THz)', 'Mode frequency', 'Hz'])
        freq_sim  = '[{:.3f} ... {:.3f}]'.format(self._dic['freq_sim'][1]*1e-12,  self._dic['freq_sim'][-2]*1e-12)
        table.add_row(['freq_sim', freq_sim + ' (THz)', 'Mode frequency to match the sim domain', 'Hz'])

        table.add_row(['Dint', '[array for each pump]', 'Integrated Dispersion at the pumps from file', 'Hz'])
        table.add_row(['Dint_sim', '[array for each pump]', 'Integrated Dispersion at the pumps from fit', 'Hz'])
        Dint0 = '[{:.3f} ... {:.3f}]'.format(self._dic['Dint0'][1]*1e-9, self._dic['Dint0'][-2]*1e-9)
        table.add_row(['Dint0', Dint0 + ' (GHz)', 'Dint at the center of sim domain', 'Hz'])
        D1 = f'{self._dic["D1"]*1e-12 :.3f}'
        table.add_row(['D1', D1, 'angular repetition rate', '(x1e12  THz)'])
        FSR = ['{:.3f}'.format(ii*1e-9) for ii in self._dic['FSR']]
        FSR = '[' + ', '.join(FSR)+ ']'
        table.add_row(['FSR', FSR + ' (GHz)', 'Free Spectra Range at the pumps', 'Hz'])
        FSR_center = '{:.3f}'.format(self._dic['FSR_center']*1e-9)
        table.add_row(['FSR_center', FSR_center + ' (GHz)', 'Free Spectra Range at the center of the domain', 'Hz'])
        table.vrules = False
        table.header = True
        table.align = "l"
        table.padding_width = 2
        table.padding_height = 25
        str_table = table.get_string()
        if self._do_tr:
            for ii in greek.items():
                str_table = str_table.replace(ii[0], ii[1] )
        return str_table

    def __repr__(self):
        if self._which == 'res':
            return self.reprRes()
        elif self._which == 'sim':
            return self.reprSim()
        elif self._which == 'sol':
            return self.reprSol()
        elif self._which == 'disp':
            return self.reprDisp()
