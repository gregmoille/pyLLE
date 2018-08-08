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


It is important to note the format of the dispersion file *TestDispersion.txt*. It be format such as each line represent a resonance, with first the azimuthal mode order then the frequency of resonance separated by a comma ',' 


|

Here, the dict keys can be defined with Greek letters for a nicer script (in my opinion) or with Latin letter such that the translator dictionary is defined by:

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

.. image:: https://github.com/gregmoille/pyLLE/blob/master/example/images/Analyze.png
    :width: 400px
    :align: center
    :alt: alternate text

The green dot corresponds to the raw data provided by the *TestDispersion.txt* file, the blue solid line is the spline fit of the integrated dispersion and the orange dashed line is the integrated dispersion that will be used in the LLE solver

|

To start the simulation, first we need to setup an hdf5 file which makes the bridge between python and julia 

.. code-block:: python
    
    solver.Setup()

Here, two attributes has been created *self.sim* and *self.res*  which are dictionaries where one can easily retrieve the different parameters of the simulation

|

For the sake of simplicity and to retrieve the difference parameter, all keys of the sim and res dictionary are translated with the previously described translator dictionary and cannot be access anymore with the Greek alphabet. Hence, in an ipython console, one can retrieve the parameter with, for example 

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

Here a attribute called *.sol* is created, corresponding to a dictionary where all the different data can be easily retrieved if needed.

|

We can finally start to plot the result of the simulation. One can start with a complete overview of the simulation, were a spectral and a temporal map Vs the LLE step is displayed in addition to the comb power Vs the LLE step

.. code-block:: python
    
    solver.PlotCombPower()

.. image:: https://github.com/gregmoille/pyLLE/blob/master/example/images/CombResults.png
    :width: 400px
    :align: center
    :alt: alternate text

From there, we can find the step of the LLE where we want to see the spectra:

.. code-block:: python
    
    ind = 600
    solver.PlotCombSpectra(ind)

.. image:: https://github.com/gregmoille/pyLLE/blob/master/example/images/CombSpectra.png
    :width: 400px
    :align: center
    :alt: alternate text

The temporal profile of the soliton can also be retrieve using 

.. code-block:: python
    
    ind = 600
    solver.PlotSolitonTime(600)

.. image:: https://github.com/gregmoille/pyLLE/blob/master/example/images/SolitonTime.png
    :width: 400px
    :align: center
    :alt: alternate text

|

One can also solver quickly the LLE through a steady state method finding the root of the LLE

.. code-block:: python
    
    sim['δω'] =  -10e9,
    solver.SolveSteadySteate()

.. image:: https://github.com/gregmoille/pyLLE/blob/master/example/images/SteadyState.png
    :width: 400px
    :align: center
    :alt: alternate text

|

Finally, an easy way to solve the simulation has been implemented in the *self.SaveResults* method where the whole class is pickled, hence can be easily retrieved following the example below:

.. code-block:: python
    
    solver.SaveResults('TestSimulation.pkl', path = 'Results/')
    import pickle as pkl
    old_solver = pkl.load(open('Results/TestSimulation.pkl','br'))


|
|

All the different method to plot the field or the spectra return the corresponding values, and one can easily access it through the *self.sol* attribute (dictionary), The complete script is: 

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

    # --  Setup the Solver --
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

