import matplotlib.pyplot as plt
import numpy as np
import pyLLE

plt.close('all')


res = {'R': 23e-6,
       'Qi': 1e6,
       'Qc': 1e6,
       'γ': 2}

sim = {'Pin': 150e-3,
       'Tscan': 2e5,
       'δω_stop': "None",
       'f_pmp': 191e12,
       'δω_init': 2e9*2*np.pi, 
       'δω': -15e9,
       'δω_end': -8e9*2*np.pi, 
       'μ_sim': [-74,170],
       'μ_fit': [-71, 180],
       'dispfile': 'TestDispersion.txt'
        }

# --  Setup thte Solver --
solver = pyLLE.LLEsovler(sim=sim,
                       res=res,
                       debug=True)
solver.Analyze(plot=True,
               plottype='all')
solver.Setup()

# --  Solver the Steady State LLE --
Ering, Ewg, f, ax = solver.SolveSteadySteate()

# # --  Solver the Temporal LLE --
# solver.SolveTemporal()
# solver.RetrieveData()
# solver.PlotCombPower()
# freq, Sout, Sring, fS, axS = solver.PlotCombSpectra(600)
# t, U, ft, axt = solver.PlotSolitonTime(600)