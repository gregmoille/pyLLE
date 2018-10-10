import numpy as np
import scipy.interpolate as itp
import scipy.io as io
import matplotlib as mpl
import matplotlib.pyplot as plt
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
import plotly.graph_objs as go

# backend = mpl.get_backend()
# print(backend)

class AnalyzeDisp(object):
    '''
    Calls to analyze the dispersion of a simulated resonator
    Initialization input. Everything is in SI ([]=facultative):
    **kwargs
        - f_center <float>: pump frequency
        - file <str>: .txt file to load 
        - R <float>: radius of the resonator
        - rM_fit <list>.: lower and upper bonds of mode to fit the dispersion
        - rM_sim <list>.: lower and upper bonds of the modes to extrapolate for the simulation
        - f <obj>: matplotlib figure handle
        - ax <obj>: matplotlib axes handle
        - label <list>: list of the string for each plot to be labelled
        - plottype <str>: define the type of plot 'all' [defaults], 'sim', 'fit', 'fem', 'ind'

    '''

    def __init__(self, **kwargs):
        self.c = 299792458
        self.dispfile = kwargs.get('file', None)
        self.R = kwargs.get('R', 23e-6)
        # self.n0 = kwargs.get('n0', 2)
        self.f_center = kwargs.get('f_center', 1550e-9)
        self.rM_fit = kwargs.get('rM_fit', [])
        self.rM_sim = kwargs.get('rM_sim', [])
        self.debug = kwargs.get('debug', False)
        self.zero_lines = kwargs.get('zero_lines', True)
        self.f = kwargs.get('f', None)
        self.ax = kwargs.get('ax', None)
        self.label = kwargs.get('label', None)
        self.plottype = kwargs.get('plottype', 'all')
        self._logger = kwargs.get('logger', None)
        self.pyType = kwargs.get('pyType', 'normal')
        assert type(self.dispfile) is str, 'Please provide a file'
        assert len(self.rM_fit) > 1, 'Please the modes to fit'

    def GetDint(self):
        '''
        Retrieve the dispersion of a resonator based on the frequency of 
        resonance and azimuthal mode order. The data are fit using a cubic spline method

        **Output**

            - self.PrM_fit: scipy.interpolate object which fitted the data
            - self.Dint_fit: fitted integrated dispersion for the simulated mode 
            - self.neff_pmp: effective index at the pump frequency
            - self.ng_pmp: group index at the pump frequency
        '''

        # -- Get Parameters --
        FSR0 = 1e12*23e-6/self.R
        R = self.R
        c0 = self.c
        rM_fit = self.rM_fit
        rM_sim = self.rM_sim
    
        lines = open(self.dispfile,'r').readlines()
        rM = np.array([float(ll.strip().split(',')[0])  for ll in lines])
        rf = np.array([float(ll.strip().split(',')[1])  for ll in lines])
        
        # -- Find the pumped mode --
        pump_index = np.where(np.abs(rf-self.f_center) < 0.5*FSR0)[0]
        assert pump_index.size == 1, 'Wavelength not found!'
        pump_index = pump_index[0]

        # -- Compute the effective/group index --
        L = 2*np.pi*R
        df = np.gradient(rf)
        self.neff = rM * c0/(rf*2*np.pi*R)
        self.neff_pmp = self.neff[pump_index]
        self.ng = c0/(df*L)
        self.ng_pmp = self.ng[pump_index]

        # Compute the dispersion
        d1_vg = np.gradient(self.ng)/c0
        D = -(rf**2/c0) * d1_vg/df
        self.D = D

        # -- Compute Dint --
        dm = np.array([-2, -1, 0, 1, 2])
        drf = rf - rf[pump_index]
        # ipdb.set_trace()
        Dfit = np.polyfit(dm, drf[pump_index+dm], 2)
        FSR_p = Dfit[1]
        D2 = Dfit[0]
        D1 = FSR_p * 2*np.pi
        mu = rM - rM[pump_index]
        w = 2*np.pi*rf
        Dint = w - (w[pump_index] + D1 * mu)

        # dm = np.arange(-2, 2.1, dtype=int)
        # drf = rf - rf[pump_index]
        # Dfit = np.polyfit(dm, drf[pump_index+dm], 2)
        # FSR_p = Dfit[1]
        β2 = -self.ng_pmp/self.c*(2*Dfit[0])/Dfit[1]**2/2/np.pi
        dφ = -(drf-(rM-rM[pump_index])*FSR_p)/FSR_p*2*np.pi

        ind_M = np.arange(
            pump_index+rM_fit[0], pump_index+rM_fit[1]+1, dtype=int)
        ind_fit = np.arange(rM_fit[0], rM_fit[1]+1, dtype=int)
        ind_sim = np.arange(rM_sim[0], rM_sim[1]+1, dtype=int)
        assert ind_M[0] > 0, 'Left range for mode order not correct'
        assert ind_M[-1] < rM.size, 'Right range for mode order not correct'

        M2fit = rM[ind_M] - rM[pump_index]
        Dint2fit = Dint[ind_M]
        dφ2fit = dφ[ind_M]
        PrM_fit = itp.splrep(M2fit, Dint2fit)
        Prφ_fit = itp.splrep(M2fit, dφ2fit)
        Dint_fit1 = itp.splev(ind_fit, PrM_fit)
        Dint_fit2 = itp.splev(ind_sim, PrM_fit)
        dφ_fitted = itp.splev(ind_sim, Prφ_fit)
        # ipdb.set_trace()

        tR = L*self.ng_pmp/c0
        # ipdb.set_trace()
        dw = np.arange(rM_sim[0], rM_sim[-1]+1)*2*np.pi/tR
        w0 = rf[pump_index] * 2*np.pi
        rf_fit = (w0+dw)/(2*np.pi)

        # ipdb.set_trace()
        if self.debug:
            if not self.pyType == 'jupyter':
                f = self.f
                ax = self.ax
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

                if self.zero_lines:
                    ax.plot([rf[0]*1e-12, rf[-1]*1e-12], [0, 0], 'k--')
                if self.plottype == 'all':
                    # ipdb.set_trace()
                    
                    
                    ax.plot(rf_fit*1e-12, Dint_fit2*1e-9/(2*np.pi), '--', lw=3, label = 'Sim')
                    ax.plot(rf[ind_M]*1e-12, Dint_fit1*1e-9/(2*np.pi), lw=1, label = 'Fit')
                    ax.plot(rf*1e-12, Dint*1e-9/(2*np.pi), '.', ms=4, alpha = 0.5, label = 'Raw')
                    
                    ax.legend()
                if self.plottype.lower() == 'sim':
                    ax.plot(rf_fit*1e-12, Dint_fit2*1e-9 /
                            (2*np.pi), label=self.label)
                    ax.legend()
                if self.plottype.lower() == 'fit':
                    ax.plot(rf[ind_M]*1e-12, Dint_fit1*1e-9/(2*np.pi), lw=2)
                    ax.plot(rf*1e-12, Dint*1e-9/(2*np.pi), '.', ms=4, alpha = 0.5)
                if self.plottype.lower() == 'fem':
                    ax.plot(rf*1e-12, Dint*1e-9/(2*np.pi), '-', label=self.label)
                    ax.legend()
                if self.plottype.lower() == 'ind':
                    ax.plot(rf[ind_M]*1e-12, Dint[ind_M]*1e-9 /
                            (2*np.pi), '-', label=self.label)
                    ax.legend()
                ax.set_xlabel('Frequency (THz)')
                ax.set_ylabel('Dint (GHz)')
                f.show()
                plt.pause(0.25)
            else:
                init_notebook_mode(connected=True)
                trace0 = go.Scatter(
                            x = rf_fit*1e-12,
                            y = Dint_fit2*1e-9/(2*np.pi),
                            mode = 'lines',
                            name = 'LLE simulation'
                        )
                trace1 = go.Scatter(
                    x = rf[ind_M]*1e-12,
                    y = Dint_fit1*1e-9/(2*np.pi),
                    line = dict(
                        width = 2,
                        dash = 'dash'),
                    name = 'Fit'
                )
                trace2 = go.Scatter(
                    x = rf*1e-12,
                    y = Dint*1e-9/(2*np.pi),
                    mode = 'markers',
                    name = 'FEM Simulation'
                )
                if self.plottype.lower() == 'all':
                    data = [trace2, trace0, trace1]
                if self.plottype.lower() == 'sim':
                    data = [trace0]
                if self.plottype.lower() == 'fit':
                    data = [trace2, trace1]
                if self.plottype.lower() == 'fem':
                    data = [trace2]

                layout = dict(xaxis = dict(title = 'Frequency (THz)'),
                  yaxis = dict(title = 'D<sub>int</sub> (GHz)'),
                  )
                fig = go.Figure(data=data, layout=layout)
                iplot(fig)
        else:
            print('coucou')


            # f, ax = plt.subplots()
            # ax.plot([ind_sim[0], ind_sim[-1]], [0, 0], 'k--')

        # f.show()
        Info = '-- Dispersion Analysis --\n'
        Info += '\tPump index: {}\n'.format(pump_index)
        Info += '\tCenter Pump: {:.3f} THz\n'.format(rf[pump_index]*1e-12)
        Info += '\tFSR: {:.2f} GHz\n'.format(FSR_p*1e-9)
        Info += '\tD2: {:.2f} MHz\n'.format(D2*1e-6)
        print(Info)

        if self._logger:
            Info += '\t mu_fit: [{:.0f}, {:.0f}]\n'.format(self.rM_fit[0],
                                                           self.rM_fit[1])
            Info += '\t mu_sim: [{:.0f}, {:.0f}]\n'.format(self.rM_sim[0],
                                                           self.rM_sim[1])

            self._logger.info('AnalyzeDisp', Info)

        self.PrM_fit = PrM_fit
        self.β2 = β2
        self.rf = rf[ind_M]
        self.rM = rM[ind_M]
        self.mu = mu[ind_M]
        self.D = self.D[ind_M]
        self.ng = self.ng[ind_M]
        self.neff = self.neff[ind_M]
        self.Dint = Dint[ind_M]
        self.Dint_fit = Dint_fit2
        self.ind_M = ind_M
        self.pump_index = pump_index
        self.dφ = dφ_fitted
        try:
            self.f = f
            self.ax = ax
        except:
            self.f = None
            self.ax = None
        return PrM_fit, Dint_fit2, self.neff_pmp, self.ng_pmp, self.f, self.ax


# if __name__ == '__main__':
#     plt.close('all')
#     analyze = AnalyzeDisp(file='TestDispersion.mat',
#                           lbd_center = 1500e-9,
#                           rM_fit=[-50, 160],
#                           debug=True)
#     analyze.GetBeta()
