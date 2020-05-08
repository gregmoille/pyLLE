import matplotlib.pyplot as plt
import numpy as np
import os
import time
import pyLLE

plt.close('all')
try:
    os.remove('LLE.log')
except:
    pass

res = {'R': 23e-6,
       'Qi': 1e6,
       'Qc': 1e6,
       'γ': 1.55,
       'dispfile': 'TestDispersion.csv'}

sim = {'Pin': 100e-3,
       'Tscan': 1e6,
       'f_pmp': 191e12,
       'δω_init': 2e9*2*np.pi,
       'δω_end': -8e9*2*np.pi,
       'μ_sim': [-74,170],
       'μ_fit': [-71, 180],
        }

# --  Setup thte Solver --
solver = pyLLE.LLEsolver(sim=sim,
                         res=res)
f = solver.Analyze(plot=True,
               plottype='all')
solver.Setup()
time.sleep(1)
# --  Solver the Temporal LLE --
solver.SolveTemporal()
time.sleep(1)
solver.RetrieveData()
solver.PlotCombPower()
ind = 450
figSpectra = solver.PlotCombSpectra(ind)
figTime = solver.PlotSolitonTime(ind)
