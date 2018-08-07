# pyLLE ![NIST logo](images/NISTlogo32x32.jpg)

![](https://readthedocs.org/projects/pylle/badge/?version=latest) 
[![](https://img.shields.io/github/license/mashape/apistatus.svg)](licence.txt)

pyLLE is a tool to solve the Lugiato Lefever Equations (LLE)<sup>[1](#ref1)</sup><sup>,</sup><sup>[2](#ref2)</sup><sup>,</sup><sup>[3](#ref3)</sup>in a fast and easy way. Thanks to an user-friendly front-end in python and a efficient back end in Julia, solving this problem becomes easy and fast. 

For a complete documentation of the package, please visit the [readthedocs page](http://pylle.readthedocs.io/en/latest/index.html)

## Instalation

As pyLLE relies on a Julia v0.6.4 back-end, please prior to install this package be sure that Julia is installed on your machine or visit the julia [package downloader page](https://julialang.org/downloads/) to install it. Please, keep julia in the default directory during the installation (i.e. ~\AppData\Local\Julia-0.6.4\). If not, please go to the manual installation. Once Julia installed, the different packages needed to run pyLLE, either python or julia related, will be automatically downloaded and installed. Just a heads up, the installation of the package can vary in time, especially because of Julia that might rebuilds the cache.

For a automatic install, just

```bash
pip install pyLLE
```

For a manual install, download the .zip of the repository or clone it and install with the setup.py script 

```bash
git clone https://github.com/gregmoille/pyLLE.git
cd pyLLE
python setup.py install
```

If the julia location is custom, please before installing change in the setup.py, line 18 to the correct location, as in pyLLE/llesolver.py line 430 to point to the correct location. Thanks

## Example

First download the .mat file available in the example folder of this repository. 


Import the pyLLE module, and the different utility modules

```python
import pyLLE
```

Define the resonator properties. Here we will simulate the ring resonator made of Si3N4 from Li et al<sup>[4](#ref4)</sup>. Hence the _res_ dictionary should look like

```python
res = {'R': 23e-6, # ring radius
       'Qi': 1e6,  # intrinsic quality factor
       'Qc': 1e6,  # coupling quality factor
       'γ': 2, # non-linear coefficient
       }
```

Now, we should define the simulation parameters through the _sim_ dictionary


```python
sim = {'Pin': 100e-3, #input power in W
       'Tscan': 0.6e5, # Total time for the simulation in unit of round trip
       'δω_stop': "None", # if we want to stop the detuning at a given frequency
       'f_pmp': 191e12, # frequency of the pump in Hz
       'δω_init': -4, # start frequency of detuning ramp in Hz
       'δω_end': 10, # stop frequency of detuning ramp in Hz
       'μ_sim': [-70,170], # limit of the mode on the left and right side of the pump to simulate
       'μ_fit': [-60, 160], # limit of the mode on the left and right side of the pump to fit the dispersion with
       'dispfile': 'ExampleQing.mat' # name of the dispersion file
        }
```


Here, the dict keys can be defined with greek letters for a nicer script (in my opinion) or with Latin letter such that the translator dictionary is defined by:

```python
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
```
Hence, instead of providing *δω_init* in the sim dictionary, one could provide *domega_init*.

We can now setup the pyLLE class: 

```python 
solver = pyLLE.LLEsovler(sim=sim,
                       res=res,
                       debug=True)
```

The debug input allows the script to generate a log file in the working directory with useful information on the simulation. The authors highly encourage to keep this key to True, unless some loops are run which could create an issue with the read/write access to the file. 


To analyze the dispersion just the *Analyze* method with the correct parameters listed in the [docs](http://pylle.readthedocs.io/en/latest/source/pyLLE.html)

```python
solver.Analyze(plot=True,
               plottype='all')
```

To start the simulation, first we need to setup an hdf5 file which makes the bridge between python and julia 

```python
solver.Setup()
```
For the sake of simplicity and to retrieve the difference parameter, all keys of the sim and res dictionary are translated with the previously described translator dictionary and cannot be access anymore with the greek alphabet. Hence, in an ipython console, one can retrieve the parameter with, for example 

```python
IN [1]: solver.sim['mu_sim']
OUT [1]: [-70,170]
```


Then we can start the simulation 

```python 
solver.Solve()
```

To retrieve the data computed by julia, we call the *RetrieveData* method

```python
solver.RetrieveData()
```

We can finally start to plot the result of the simulation. One can start with a complete overview of the simulation, were a spectral and a temporal map Vs the LLE step is displayed in addition to the comb power Vs the LLE step

```python
solver.PlotCombPower()
```

From there, we can find the step of the LLE where we want to see the spectra:

```python
ind = 600
solver.PlotCombSpectra(ind)
```


One can also solver quickly the LLE through a steady state method finding the root of the LLE

```python
sim['δω'] =  -10e9,
Ering, Ewg, f, ax = solver.SolveSteadySteate()
```



The complete script is: 

```python
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
       'dispfile': 'TestDispersion.mat'
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
```

## How to Cite Us?

Soon you will be. For the moment, please provide the name of the package, the authors (Gregory Moille, Qing Li, Xiyuan Lu and Kartik Srinivasan) as a full url to the repository


## References

<a name="ref1">1</a>: Luigi A. Lugiato and René Lefever. "Spatial dissipative structures in passive optical systems." Physical review letters 58, no. 21 (1987): 2209.

<a name="ref1">2</a>: Yanne K. Chembo and Curtis R. Menyuk. "Spatiotemporal Lugiato-Lefever formalism for Kerr-comb generation in whispering-gallery-mode resonators." Physical Review A 87, no. 5 (2013): 053852.

<a name="ref1">3</a>: Stéphane Coen, Hamish G. Randle, Thibaut Sylvestre, and Miro Erkintalo. "Modeling of octave-spanning Kerr frequency combs using a generalized mean-field Lugiato–Lefever model." Optics letters 38, no. 1 (2013): 37-39.

<a name="ref1">4</a>: Qing Li, Travis C. Briles, Daron A. Westly, Tara E. Drake, Jordan R. Stone, B. Robert Ilic, Scott A. Diddams, Scott B. Papp, and Kartik Srinivasan. "Stably accessing octave-spanning microresonator frequency combs in the soliton regime." Optica 4, no. 2 (2017): 193-203.