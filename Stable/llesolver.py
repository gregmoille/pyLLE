import numpy as np
from .analyzedisp import AnalyzeDisp
import scipy.interpolate as itp
import matplotlib.pyplot as plt
import time
import sys
import scipy
import subprocess as sub
import os
from copy import copy
import shutil
import tempfile
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)
    import h5py

import ipdb
# path_juliaScript = __file__.split['/']

path_juliaScript = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/'
tmp_dir = tempfile.gettempdir() + '/'
# print(tmp_dir)
# print('-'*30)
# print(path_juliaScript)


class LLEsovler(object):
    _c0 = 299792458
    _hbar = 6.634e-34/(2*np.pi)

    def __init__(self, **kwargs):
        self.res = kwargs.get('res', {})
        self.sim = kwargs.get('sim', {})
        self.sim_norm = kwargs.get('sim_norm', {})
        self._debug = kwargs.get('debug', False)

        # find all the needed parameters
        assert 'Qi' in self.res.keys(), 'Please provide Qi'
        assert 'Qc' in self.res.keys(), 'Please provide Qc'
        assert 'R' in self.res.keys(), 'Please provide R'
        # assert 'V' in self.res.keys(), 'Please provide V'
        # assert 'Pin' in self.sim.keys(), 'Please provide Pin'
        assert 'Tscan' in self.sim.keys(), 'Please provide Tscan'
        # assert 'dphi_init' in self.sim.keys(), 'Please provide dphi_init'
        # assert 'dphi_end' in self.sim.keys(), 'Please provide dphi_end'
        # assert 'lbd_pmp' in self.sim.keys(), 'Please provide lbd_pmp'
        assert 'mu_sim' in self.sim.keys(), 'Please provide mu_sim'
        assert 'mu_fit' in self.sim.keys(), 'Please provide mu_fit'
        assert 'dispfile' in self.sim.keys(), 'Please provide dispfile'

    def Analyze(self, plot=False, f=None, ax=None, label=None, plottype='all', zero_lines = True):
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

        do_plot = self._debug or plot
        self._analyze = AnalyzeDisp(file=self.sim['dispfile'],
                                    f_center=self.sim['f_pmp'],
                                    rM_fit=self.sim['mu_fit'],
                                    rM_sim=self.sim['mu_sim'],
                                    R=self.res['R'],
                                    debug=do_plot,
                                    f=f,
                                    ax=ax,
                                    label=label,
                                    plottype=plottype,
                                    zero_lines = zero_lines)
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

    def Setup(self):
        # -- Make the hdf5 file --
        # ------------------------------------------------------------
        try:
            os.remove(tmp_dir + 'ParamLLEJulia.h5')
            os.remove(tmp_dir + 'ResultsJulia.h5')
        except:
            pass

        # -- Normalize parameters --
        # ------------------------------------------------------------
        w0 = 2*np.pi * self.sim['f_pmp']
        hbar = self._hbar
        µ_sim = self.sim['mu_sim']
        c0 = self._c0
        ng = self.res['ng']
        keys = self.sim_norm.keys()
        µ = np.arange(µ_sim[0], µ_sim[-1]+1)
        pmp_ind = np.where(µ==0)[0][0]
        # -- Normalized coupling losses --
        if not('dw_ext' in keys):
            Qc = self.res['Qc']
            Q0 = self.res['Qi']
            self.sim_norm['dw_ext'] = w0 * (1/Qc)

        # -- Normalized total losses --
        if not('dw_tot' in keys):
            Qc = self.res['Qc']
            Q0 = self.res['Qi']
            self.sim_norm['dw_tot'] = w0 * \
                (1/self.res['Qc'] + 1/self.res['Qi'])

        # -- Normalized Pump Power --
        if not('F2' in keys):
            n2 = self.res['n2']
            V = self.res['V']
            Pin = self.sim['Pin']
            dw_ext = self.sim_norm['dw_ext']
            dw_tot = self.sim_norm['dw_tot']
            g0 = n2 * c0 * hbar * w0**2/(ng**2 * V)
            fact_normE = np.sqrt(8*g0*dw_ext/(dw_tot**3 * hbar * w0))
            F2 = Pin * fact_normE**2
            self.sim_norm['g0'] = g0
            self.sim_norm['fact_normE'] = fact_normE
            self.sim_norm['F2'] = F2
        
            
        else:
            if not ('fact_normE' in self.sim_norm):
                n2 = self.res['n2']
                V = self.res['V']
                Pin = self.sim['Pin']

                dw_ext = self.sim_norm['dw_ext']
                dw_tot = self.sim_norm['dw_tot']
                g0 = n2 * c0 * hbar * w0**2/(ng**2 * V)
                fact_normE = np.sqrt(8*g0*dw_ext/(dw_tot**3 * hbar * w0))
                self.sim_norm['fact_normE'] = fact_normE
            if not('g0' in keys):
                n2 = self.res['n2']
                V = self.res['V']
                g0 = n2 * c0 * hbar * w0**2/(ng**2 * V)
                self.sim_norm['g0'] = g0

        # -- Normalized Detunning --
        if not('alpha_init' in keys):
            self.sim_norm['alpha_init'] = -2*self.sim['dw_init']/dw_tot
        if not('alpha_end' in keys):
            self.sim_norm['alpha_end'] = -2*self.sim['dw_end']/dw_tot
        if not('alpha_stop' in keys):
            if not('dw_stop' in self.sim.keys()):
                self.sim_norm['alpha_stop'] = "None"
            else:
                self.sim_norm['alpha_stop'] = self.sim['dw_stop']
        
        # -- Coupling Dispersion -- 
        if not('domega' in keys):
            if not ('dQc' in self.res.keys()):
                domega = lambda µ: 1 + 0*µ
            else:
                domega = lambda µ: 1/self.res.dQc(µ)
        else:
            domega = self.sim_norm['domega']
        
        domega_val = domega(µ)
        self.sim_norm['domega_val'] = domega_val/domega_val[pmp_ind]
        # -- Misc --
        L =  2*np.pi*self.res['R']
        self.sim_norm['Tscan'] = self.sim['Tscan']
        self.sim_norm['mu_sim'] = self.sim['mu_sim']
        tR = L *self.res['ng']/self._c0
        self.sim_norm['tR'] = tR
        self.sim_norm['Dint'] =  self.sim['Dint']
        self.sim_norm['t_end'] = self.sim['Tscan']*tR*self.sim_norm['dw_tot']/2
        w = w0 + µ *2*np.pi/tR
        self.sim_norm['w'] = w
        self.sim_norm['mu'] = µ
        # -- Print the Normalized Parameters --
        print('Normalized Parameters:')
        for k, it in self.sim_norm.items():
            if not (k == 'w' or k =='Dint' or k == 'domega_val' or k == 'mu'):
                if not type(it) is list:
                    try:
                        print('\t{} = {:.2f}'.format(k, it))
                    except:
                        print('\t{} = {}'.format(k, it))
                else:
                    print('\t{} = {}'.format(k, it))

        # -- create h5file --
        h5f = h5py.File(tmp_dir + 'ParamLLEJulia.h5', 'w')
        h5f.create_group('sim_norm')
        cnt = 0
        for key, it in self.sim_norm.items():
            if not key == 'domega':
                if type(it) is str:
                    it = np.string_(it)
                # print('sim_norm/{}'.format(key))
                h5f.create_dataset('sim_norm/{}'.format(key), data=[it])
        cnt += 1

        h5f.close()

    def Solve(self):
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

    def RetrieveData(self):
        drct = tmp_dir
        S = h5py.File(tmp_dir + 'ResultsJulia.h5', 'r')
        sol = {}
        keys = ['u_probe','Em_probe', 'Ewg', 'comb_power']
        for k in keys:
            rl = 'Results/{}Real'.format(k)
            im = 'Results/{}Imag'.format(k)
            sol[k] = S[rl][:] + 1j*S[im][:]
        sol['detuning'] = S["Results/detuningReal"][:]
        S.close()
        os.remove(tmp_dir + 'ParamLLEJulia.h5')
        os.remove(tmp_dir + 'ResultsJulia.h5')
        sol['freq'] = self.sim_norm['w']/(2*np.pi)

        # Normalize to SI
        g0 = self.sim_norm['g0'] 
        dw_tot = self.sim_norm['dw_tot'] 
        # normE = np.sqrt(2*g0/dw_tot)
        normE = self.sim_norm['fact_normE']
        beta = self.sim_norm['Dint']/dw_tot
        Nprobe = sol['Em_probe'].shape[1]
        t_end = self.sim_norm['t_end']
        mu_sim = self.sim['mu_sim']
        tau = np.linspace(0,t_end, Nprobe)

        mu = np.arange(mu_sim[0], mu_sim[-1]+1)


        cnt = 0
        # for tt in tau: 
        #     phi = np.exp(-1j*beta*tt)
        #     Epb = sol['Em_probe'][:,cnt]
        #     Ewg = sol['Ewg'][:,cnt]
        #     sol['Em_probe'][:,cnt] = phi*np.conj(Epb)/normE
        #     sol['Ewg'][:,cnt] = phi*np.conj(Ewg)/normE
        #     th_cnt = 0
        sol['theta'] = np.linspace(-np.pi,np.pi, sol['u_probe'].shape[0])
        sol['u_probe'] = np.conj(sol['u_probe'])/normE
        self.sol = sol

    def PlotCombPower(self):
        cmap = plt.get_cmap("tab10")
        freq = self.sol['freq']*1e-12
        Epb = self.sol['Em_probe']
        Epb = 10*np.log10(np.abs(Epb)**2)
        Epb = Epb - Epb.max()
        # Epb[Epb<=0] = 1e-100
        f, ax = plt.subplots(3, 1, sharex=True)
        ax[0].pcolormesh(np.arange(0, 1000), freq,Epb,
                         rasterized=True)
        theta = self.sol['theta']
        ax[1].pcolormesh(np.arange(0, 1000), theta,
                         (np.abs(self.sol['u_probe'])**2),
                         rasterized=True, vmin=0)
        ax[2].plot(np.arange(0, 1000), self.sol['comb_power'] /
                   self.sol['comb_power'].max())
        ax = np.append(ax, ax[2].twinx())
        ax[3].plot(np.arange(0, 1000),self.sol['detuning'],
                    c = cmap.colors[1])
        ax[0].set_ylabel('Frequency (THz)')
        ax[1].set_ylabel('Angle (x π)')
        ax[2].set_xlabel('LLE Step (sub-sampled)')
        ax[2].set_ylabel('Norm. Comb Pwr')
        ax[3].set_ylabel('Norm. Detuning')
        ax[0].set_xlim([0,1000])
        f.show()

        return f, ax

    def PlotCombSpectra(self, ind, f=None, ax=None, label=None, pwr='both'):
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

        Sring = 10*np.log10(1e3*np.abs(self.sol['Em_probe'][:, ind])**2)
        Sout = 10*np.log10(1e3*np.abs(self.sol['Ewg'][:, ind])**2)

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

        return f, ax, freq, Sout, Sring

    def PlotSolitonTime(self, ind):
        pass

    def SaveResults(self, fname, path='./'):
        pass

# if __name__ == '__main__':
#     plt.close('all')
#     sim = {'Pin': 200e-3,
#            'Tscan': 1e5,
#            'dphi_init': -10,
#            'dphi_end': [5, 1],
#            'lbd_pmp': 1069.5e-9,
#            'mu_sim': [-150, 170],
#            'mu_fit': [-100, 100],
#            'dispfile': 'Data/AcesRingBoulder.mat'}
#     res = {'R': 23e-6,
#            'Qi': 1e6,
#            'Qc': 1e6,
#            'gamma': 2, }

#     solver = SolveLLE(sim=sim, res=res, debug=False)
#     solver.Analyze(plot=True)
#     # plt.pause(0.15)
#     solver.Setup()
#     solver.Solve()
#     solver.RetrieveData()
#     solver.PlotCombPower()
