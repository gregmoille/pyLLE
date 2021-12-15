import numpy as np
import scipy.interpolate as itp
import scipy.io as io
from scipy import constants as cts
import matplotlib as mpl
import matplotlib.pyplot as plt
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
import plotly.graph_objs as go
import ipdb
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
        cts.c = 299792458
        self.dispfile = kwargs.get('file', None)
        self.R = kwargs.get('R', 23e-6)
        self.f_center = kwargs.get('f_center', None)
        self.fpmp = kwargs.get('fpmp', None)
        self.rM_fit = kwargs.get('rM_fit', [])
        self.rM_sim = kwargs.get('rM_sim', [])
        self.debug = kwargs.get('debug', False)
        self.plottype = kwargs.get('plottype', 'all')
        self._logger = kwargs.get('logger', None)
        self.pyType = kwargs.get('pyType', 'normal')
        self.fig = kwargs.get('fig', 'normal')
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


        def getFreqModes(f_center= None):
            FSR0 = 1e12*23e-6/self.R
            lines = open(self.dispfile,'r').readlines()
            self.rM = np.array([float(ll.strip().split(',')[0])  for ll in lines])
            self.rf = np.array([float(ll.strip().split(',')[1])  for ll in lines])
            # -- Find the pumped mode --
            if not f_center:
                pmp_ind = np.where(np.abs(self.rf-self.fpmp[0]) < 0.5*FSR0)[0]
                assert pmp_ind.size == 1, 'Wavelength not found!'
                self.pmp_ind_fit = pmp_ind[0]
            else:
                pmp_ind = np.where(np.abs(self.rf-f_center) < 0.5*FSR0)[0]
                assert pmp_ind.size == 1, 'Wavelength not found!'
                return pmp_ind[0]
            return self.pmp_ind_fit

        def getGroupIndex():
            # -- Compute the effective/group index --
            L = 2*np.pi*self.R
            df = np.gradient(self.rf)
            self.neff = self.rM * cts.c/(self.rf*2*np.pi*self.R)
            self.neff_pmp = self.neff[self.pmp_ind_fit]
            self.ng = cts.c/(df*L)
            self.ng_pmp = self.ng[self.pmp_ind_fit]
            self.tR = L*self.ng_pmp/cts.c

        def getDispersion():
            # Compute the dispersion
            df = np.gradient(self.rf)
            d1_vg = np.gradient(self.ng)/cts.c
            D = -(self.rf**2/cts.c) * d1_vg/df
            self.D = D

        def getIntegratedDisp(ind_pmp):
            # -- Compute Dint --
            dm = np.array([-2, -1, 0, 1, 2])
            drf = self.rf - self.rf[self.pmp_ind_fit]
            Dfit = np.polyfit(dm, drf[ind_pmp+dm], 2)
            self.FSR = Dfit[1]
            self.D2 = Dfit[0]*2*np.pi
            D1 = self.FSR * 2*np.pi
            μ = self.rM - self.rM[ind_pmp]
            ω = 2*np.pi*self.rf
            Dint = ω - (ω[ind_pmp] + D1 * μ)
            β2 = -self.ng_pmp/cts.c*(2*Dfit[0])/Dfit[1]**2/2/np.pi
            self.μ0 = μ
            return μ, Dint, D1

        def fitDintDomain(ind0,ind_master, Dint):

            if self.rM_fit == [None, None]:
                μfit = [self.μ0[0], self.μ0[-1]]
                shift = 0
            else:
                μfit = self.rM_fit
                shift = ind0 - ind_master
            # print(shift)
            μfit  = [rm-shift for rm in μfit]
            # print(μfit)
            M = np.arange(ind0+μfit[0],
                          ind0+μfit[1]+1,
                          dtype=int)
            # print(M)
            μselect = np.arange(μfit[0],μfit[1]+1, dtype=int)
            assert M[0] >= 0, 'Left range for mode order not correct'
            assert M[-1] <= self.rM.size, 'Right range for mode order not correct'


            # self.PrM_fit = itp.splrep(μfit, Dint2fit)
            M2fit = self.rM[M] - self.rM[ind0]
            Dint2fit = Dint[M]
            fitfun = itp.splrep(M2fit, Dint2fit)
            fit_selected = itp.splev(μselect, fitfun)
            pmp_ind = np.argwhere(M2fit == 0).flatten()[0]
            return fitfun, μselect, fit_selected, pmp_ind

        def doFitSimulationDint(fitfun, ind0, ind_master):
            # if self.rM_fit == [None, None]:
            #     shift = 0
            # else:
            #     μfit = self.rM_fit
            shift = ind0 - ind_master

            μ2fit = [rm-shift for rm in self.rM_sim]
            ind_sim = np.arange(μ2fit[0], μ2fit[1]+1, dtype=int)

            Dint_fit = itp.splev(ind_sim, fitfun)

            pmp_ind = np.argwhere(ind_sim == 0).flatten()[0]

            return μ2fit, Dint_fit, pmp_ind



        # --  Get all parameters that are not relative to the pump mode --
        # ----------------------------------------------------------------------------


        self.μsim = []
        self.Dint_fit = []
        self.Dint_sim = []
        self.pmp_ind = []
        self.Dint = []
        self.D1 = []
        self.μsim = []
        self.μfit = []
        self.ind_pmp_sim = []
        self.ind_pmp_fit = []
        # --  Get all parameters that are relative to the pump mode --
        # ----------------------------------------------------------------------------
        # The data collected here are only entended for visulation purpose as the
        # integration Dint that is use in pyLLE is computed at the center of the
        # domain. For more information, please look at ComputeLLE.jl file (core
        # function in Julia of pyLLE) or chekc out Taheri et al. eq 12
        cnt_f = 0
        for ff in self.fpmp:
            if not ff == self.fpmp[0]:
                ind_pmp = getFreqModes(f_center= ff)
            else:
                ind_pmp = getFreqModes()
                ind_master = ind_pmp
                getGroupIndex()
                getDispersion()

                # will do the fit relative to the first pump
            μ_, Dint_, D1_ = getIntegratedDisp(ind_pmp)
            ff, μfit_, Dint_fit_, pmp_fit_ = fitDintDomain(ind_pmp,ind_master, Dint_)
            μsim_, Dint_sim_, pmp_sim_ = doFitSimulationDint(ff, ind_pmp, ind_master)

            self.Dint += [Dint_]
            self.Dint_fit += [Dint_fit_]
            self.Dint_sim += [Dint_sim_]
            self.D1 += [D1_]
            self.ind_pmp_sim += [pmp_sim_]
            self.ind_pmp_fit += [pmp_fit_]
            self.μsim += [μsim_]
            self.μfit += [μfit_]
            self.fpmp[cnt_f] = self.rf[ind_pmp]

            cnt_f += 1

        ff0 = self.fpmp[0]
        # print(μdomain)
        ind0 = np.sum(self.μsim[0])/2
        # ind1 = np.sum(self.μsim[1])/2
        assert ind0 == int(ind0), 'Domain not correct'
        ind_center = int(self.pmp_ind_fit + ind0)


        for ii in range(len(self.fpmp)):
            self.pmp_ind += [int(-1*np.sum(self.μsim[ii])/2)]


        f_center = self.rf[ind_center]
        ind_center = getFreqModes(f_center= f_center)
        μ_, Dint_, D1_ = getIntegratedDisp(ind_center)
        ff, μfit_, Dint_fit_, pmp_fit_ = fitDintDomain(ind_center,ind_center, Dint_)
        μsim_, Dint_sim_, pmp_sim_ = doFitSimulationDint(ff, ind_center, ind_center-ind0)

        self.fpmp += [f_center]
        self.Dint += [Dint_]
        self.Dint_fit += [Dint_fit_]
        self.Dint_sim += [Dint_sim_]
        self.D1 += [D1_]
        self.ind_pmp_sim += [pmp_sim_]
        self.ind_pmp_fit += [pmp_fit_]
        self.μsim += [μsim_]
        self.μfit += [μfit_]

        # μfit = self.μfit[ii] #+ self.ind_pmp_fit[ii]
        μsim = self.μsim[ii] #+ self.ind_pmp_sim[ii]
        self.freq_fit = self.fpmp[0] + np.arange(μsim[0], μsim[-1]+1)*self.D1[0]/(2*np.pi)

        #         dφ = -(drf-(rM-rM[self.pmp_ind_fit])*FSR_p)/FSR_p*2*np.pi

    def DisplayParam(self):

        Info = '-- Dispersion Analysis --\n'
        Info += '\tPump index: {}\n'.format(self.pmp_ind_fit)
        Info += '\tCenter Pump: {:.3f} THz\n'.format(self.rf[self.pmp_ind_fit]*1e-12)
        Info += '\tFSR: {:.2f} GHz\n'.format(self.FSR*1e-9)
        Info += '\tD2: {:.2f} MHz\n'.format(self.D2*1e-6)
        print(Info)

        if self._logger:
            Info += '\t mu_fit: [{:.0f}, {:.0f}]\n'.format(self.rM_fit[0],
                                                           self.rM_fit[1])
            Info += '\t mu_sim: [{:.0f}, {:.0f}]\n'.format(self.rM_sim[0],
                                                           self.rM_sim[1])

            self._logger.info('AnalyzeDisp', Info)

    def DisplayPlot(self, display_center = True):

        def PlotJupyter():
            init_notebook_mode(connected=True)
            trace = []
            for ii in range(len(self.fpmp)-1):
                μfit = self.μfit[ii] #+ self.ind_pmp_fit[ii]
                μsim = self.μsim[ii] #+ self.ind_pmp_sim[ii]
                dν_fit = np.arange(μfit[0], μfit[-1]+1)*self.D1[0]/(2*np.pi)
                dν_sim = np.arange(μsim[0], μsim[-1]+1)*self.D1[0]/(2*np.pi)
                ν0 = self.fpmp[ii]
                rf = self.rf
                rf_fit = (ν0+dν_fit)
                rf_sim = (ν0+dν_sim)

                trace += [go.Scatter(
                            x = rf*1e-12,
                            y = self.Dint[ii]*1e-9/(2*np.pi),
                            legendgroup = 'Pump #{:.0f}'.format(ii),
                            mode = 'markers',
                            name = 'FEM Simulation')]
                if self.plottype.lower() == 'all':

                    trace += [go.Scatter(
                                x = rf_fit*1e-12,
                                y = self.Dint_fit[ii]*1e-9/(2*np.pi),
                                line = dict(width = 2, dash = 'dash'),
                                legendgroup = 'Pump #{:.0f}'.format(ii),
                                name = 'Fit')]
                    trace += [go.Scatter(
                                x = rf_sim*1e-12,
                                y = self.Dint_sim[ii]*1e-9/(2*np.pi),
                                legendgroup = 'Pump #{:.0f}'.format(ii),
                                mode = 'lines',
                                name = 'LLE simulation')]


                if self.plottype.lower() == 'sim':

                    trace += [go.Scatter(
                                x = rf_sim*1e-12,
                                y = self.Dint_sim[ii]*1e-9/(2*np.pi),
                                legendgroup = 'Pump #{:.0f}'.format(ii),
                                mode = 'lines',
                                name = 'LLE simulation')]

                if self.plottype.lower() == 'fit':

                    trace += [go.Scatter(
                                x = rf_fit*1e-12,
                                y = self.Dint_fit[ii]*1e-9/(2*np.pi),
                                line = dict(width = 2, dash = 'dash'),
                                legendgroup = 'Pump #{:.0f}'.format(ii),
                                name = 'Fit')]
                if self.plottype.lower() == 'fem':
                    pass


            if display_center:
                ii = len(self.fpmp)-1
                μfit = self.μfit[ii] #+ self.ind_pmp_fit[ii]
                μsim = self.μsim[ii] #+ self.ind_pmp_sim[ii]
                dν_fit = np.arange(μfit[0], μfit[-1]+1)*self.D1[0]/(2*np.pi)
                dν_sim = np.arange(μsim[0], μsim[-1]+1)*self.D1[0]/(2*np.pi)
                ν0 = self.fpmp[ii]
                rf = self.rf
                rf_fit = (ν0+dν_fit)
                rf_sim = (ν0+dν_sim)

                trace += [go.Scatter(
                            x = rf*1e-12,
                            y = self.Dint[ii]*1e-9/(2*np.pi),
                            legendgroup = 'Pump #{:.0f}'.format(ii),
                            mode = 'markers',
                            name = 'FEM Simulation')]
                if self.plottype.lower() == 'all':

                    trace += [go.Scatter(
                                x = rf_fit*1e-12,
                                y = self.Dint_fit[ii]*1e-9/(2*np.pi),
                                line = dict(width = 2, dash = 'dash'),
                                legendgroup = 'Pump #{:.0f}'.format(ii),
                                name = 'Fit')]
                    trace += [go.Scatter(
                                x = rf_sim*1e-12,
                                y = self.Dint_sim[ii]*1e-9/(2*np.pi),
                                legendgroup = 'Pump #{:.0f}'.format(ii),
                                mode = 'lines',
                                name = 'LLE simulation')]


                if self.plottype.lower() == 'sim':

                    trace += [go.Scatter(
                                x = rf_sim*1e-12,
                                y = self.Dint_sim[ii]*1e-9/(2*np.pi),
                                legendgroup = 'Pump #{:.0f}'.format(ii),
                                mode = 'lines',
                                name = 'LLE simulation')]

                if self.plottype.lower() == 'fit':

                    trace += [go.Scatter(
                                x = rf_fit*1e-12,
                                y = self.Dint_fit[ii]*1e-9/(2*np.pi),
                                line = dict(width = 2, dash = 'dash'),
                                legendgroup = 'Pump #{:.0f}'.format(ii),
                                name = 'Fit')]

                trace += [go.Scatter(
                            x = [ν0*1e-12],
                            y = [0],
                            mode = 'markers',
                            marker = dict(color = 'black', size = 8),
                            name = 'center domain')]
            self.fig = go.FigureWidget(data=trace)
            self.fig.update_layout(xaxis = dict(title = 'Frequency (THz)'),
                              yaxis = dict(title = 'D<sub>int</sub> (GHz)'))

        def PlotMatplotlib():
            self.fig, ax = plt.subplots()
            for ii in range(len(self.fpmp)-1):
                μfit = self.μfit[ii] #+ self.ind_pmp_fit[ii]
                μsim = self.μsim[ii] #+ self.ind_pmp_sim[ii]
                dν_fit = np.arange(μfit[0], μfit[-1]+1)*self.D1[0]/(2*np.pi)
                dν_sim = np.arange(μsim[0], μsim[-1]+1)*self.D1[0]/(2*np.pi)
                ν0 = self.fpmp[ii]
                rf = self.rf
                rf_fit = (ν0+dν_fit)
                rf_sim = (ν0+dν_sim)

                ax.plot(rf*1e-12,
                       self.Dint[ii]*1e-9/(2*np.pi),
                       'o',ms= 3,
                       label = 'FEM Simulation')

                if self.plottype.lower() == 'all':

                    ax.plot(rf_fit*1e-12,
                           self.Dint_fit[ii]*1e-9/(2*np.pi),
                           '--',
                           label = 'Fit')
                    ax.plot(rf_sim*1e-12,
                           self.Dint_sim[ii]*1e-9/(2*np.pi),
                           label = 'LLE simulation')

                    if self.plottype.lower() == 'sim':
                        ax.plot(rf_sim*1e-12,
                               self.Dint_sim[ii]*1e-9/(2*np.pi),
                               label = 'LLE simulation')

                    if self.plottype.lower() == 'fit':
                        ax.plot(rf_fit*1e-12,
                               self.Dint_fit[ii]*1e-9/(2*np.pi),
                               label = 'Fit')

        if not self.pyType == 'jupyter':
            PlotMatplotlib()
        else:
            PlotJupyter()

        return self.fig
# if __name__ == '__main__':
#     plt.close('all')
#     analyze = AnalyzeDisp(file='TestDispersion.mat',
#                           lbd_center = 1500e-9,
#                           rM_fit=[-50, 160],
#                           debug=True)
#     analyze.GetBeta()
