import matplotlib.pyplot as plt
import pyLLE

plt.close('all')


res = {'R': 23e-6,
       'Qi': 1e6,
       'Qc': 1e6,
       'γ': 2}

sim = {'Pin': 100e-3,
       'Tscan': 0.6e5,
       'δω_stop': "None",
       'f_pmp': 191e12,
       'δω_init': -4, 
       'δω_end': 10, 
       'μ_sim': [-70,170],
       'μ_fit': [-70, 170],
       'dispfile': 'ExampleQing.mat'
        }

solver = pyLLE.LLEsovler(sim=sim,
                       res=res,
                       debug=True)
solver.Analyze(plot=True,
               plottype='all')

solver.Setup()
solver.Solve()
solver.RetrieveData()
solver.PlotCombPower()
