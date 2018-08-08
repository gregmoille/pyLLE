Quick Tutorial to simulate a ring resonator
===========================================

First download the .txt file available in the example folder of the repository. 

Import the pyLLE module, and the different utility modules

.. code-block:: python
    
    import pyLLE


Define the resonator properties. Here we will simulate the ring resonator made of Si3N4 from Li et al<sup>[4](#ref4)</sup>. Hence the _res_ dictionary should look like

.. code-block:: python

    res = {'R': 23e-6, # ring radius
           'Qi': 1e6,  # intrinsic quality factor
            'Qc': 1e6,  # coupling quality factor
            'γ': 2, # non-linear coefficient
            }


Now, we should define the simulation parameters through the _sim_ dictionary


.. code-block:: python

    sim = {'Pin': 100e-3,  #input power in W
           'Tscan': 0.6e5, # Total time for the simulation in unit of round trip
           'δω_stop': "None", # if we want to stop the detuning at a given frequency
           'f_pmp': 191e12, # frequency of the pump in Hz
           'δω_init': -4, # start frequency of detuning ramp in Hz
           'δω_end': 10, # stop frequency of detuning ramp in Hz
           'μ_sim': [-70,170], # limit of the mode on the left and right side of the pump to simulate
           'μ_fit': [-60, 160], # limit of the mode on the left and right side of the pump to fit the dispersion with
           'dispfile': 'TestDispersion.txt' # name of the dispersion file
            }


It is important to note the format of the dipersion file *TestDispersion.txt*. It be format such as each line represent an resoanace, with first the azymuthal mode order then the frequency of resonance separated by a comma ',' 


|

Here, the dict keys can be defined with greek letters for a nicer script (in my opinion) or with Latin letter such that the translator dictionary is defined by:

.. code-block:: python

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

Hence, instead of providing *δω_init* in the sim dictionary, one could provide *domega_init*.

|

We can now setup the pyLLE class: 

.. code-block:: python
    
    solver = pyLLE.LLEsovler(sim=sim,
                           res=res,
                           debug=True)

The debug input allows the script to generate a log file in the working directory with useful information on the simulation. The authors highly encourage to keep this key to True, unless some loops are run which could create an issue with the read/write access to the file. 

|

To analyze the dispersion just the *Analyze* method with the correct parameters listed in the [docs](http://pylle.readthedocs.io/en/latest/source/pyLLE.html)

.. code-block:: python
    
    solver.Analyze(plot=True,
                   plottype='all')


To start the simulation, first we need to setup an hdf5 file which makes the bridge between python and julia 

.. code-block:: python
    
    solver.Setup()


|

For the sake of simplicity and to retrieve the difference parameter, all keys of the sim and res dictionary are translated with the previously described translator dictionary and cannot be access anymore with the greek alphabet. Hence, in an ipython console, one can retrieve the parameter with, for example 

.. code-block:: python
    
    IN [1]: solver.sim['mu_sim']
    OUT [1]: [-70,170]

|

Then we can start the simulation 

.. code-block:: python
    
    solver.Solve()


To retrieve the data computed by julia, we call the *RetrieveData* method

.. code-block:: python
    
    solver.RetrieveData()


We can finally start to plot the result of the simulation. One can start with a complete overview of the simulation, were a spectral and a temporal map Vs the LLE step is displayed in addition to the comb power Vs the LLE step

.. code-block:: python
    
    solver.PlotCombPower()


From there, we can find the step of the LLE where we want to see the spectra:

.. code-block:: python
    
    ind = 600
    solver.PlotCombSpectra(ind)



One can also solver quickly the LLE through a steady state method finding the root of the LLE

.. code-block:: python
    
    sim['δω'] =  -10e9,
    Ering, Ewg, f, ax = solver.SolveSteadySteate()

|
|

The complete script is: 

.. code-block:: python

    import matplotlib.pyplot as plt
    import numpy as np
    import pyLLE

    plt.close('all')


    res = {'R': 23e-6,
           'Qi': 1e6,
           'Qc': 1e6,
           'γ': 2}

    sim = {'Pin': 100e-3,
           'Tscan': 2e5,
           'δω_stop': "None",
           'f_pmp': 191e12,
           'δω_init': 1e9*2*np.pi, 
           'δω': -10e9,
           'δω_end': -5e9*2*np.pi, 
           'μ_sim': [-70,170],
           'μ_fit': [-71, 170],
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

    # --  Solver the Temporal LLE --
    solver.SolveTemporal()
    solver.RetrieveData()
    solver.PlotCombPower()
    freq, Sout, Sring, fS, axS = solver.PlotCombSpectra(600)
    t, U, ft, axt = solver.PlotSolitonTime(600)

