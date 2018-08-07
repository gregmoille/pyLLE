import matplotlib.pyplot as plt
import numpy as np
import pyLLE

plt.close('all')


res = {'R': 23e-6,
       'Qi': 1e6,
       'Qc': 1e6,
       'γ': 2}

sim = {'Pin': 50e-3,
       'Tscan': 0.6e5,
       'δω_stop': "None",
       'f_pmp': 191e12,
       'δω_init': 2e9*2*np.pi, 
       'δω': -10e9,
       'δω_end': -5e9*2*np.pi, 
       'μ_sim': [-70,170],
       'μ_fit': [-71, 170],
       'dispfile': 'TestDispersion.mat'
        }

solver = pyLLE.LLEsovler(sim=sim,
                       res=res,
                       debug=True)
solver.Analyze(plot=True,
               plottype='all')

solver.Setup()


Ering, Ewg, f, ax = solver.SolveSteadySteate()

# solver.SolveTemporal()
# solver.RetrieveData()
# solver.PlotCombPower()
