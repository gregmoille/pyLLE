---
layout: page
title: Example
sidebar_link: true
---

## Import and Setup 

Let’s start by importing the package

```python 
import pyLLE
```

We now define the resonator parameters.

```python 
res = {'R': 23e-6, # ring radius in meter
           'Qi': 1e6,  # Intrinsic Q factor
           'Qc': 1e6,  # Coupled Q factor
           'γ': 1.55,  # Non-linear coefficient at the pump frequency
           'dispfile': 'TestDispersion.txt', # frequency and corresponding azymuthal mode simulated previously
          }
```

We now define the simulation parameters. Here we precise a linear detuning ramp of the pump from *δω_init* to *δω_end* relative to the pump mode angular frequency, mode closest to the defined pump frequency *f_pmp*. The simulation length *Tscan* is in unit of round trip, as it is more convenient in the Lugiato-Lefever formalism. It is important to notice that two parameters for the mode bandwidth have to be defined, *μ_fit* which determined the fit window of the raw data found in *dispfile*, and *μ_sim* which is the number of mode simulated in the LLE, hence could be larger than the fit mode through extrapolation

```python
import numpy as np
sim = {'Pin': 150e-3, # Input power in Q
       'Tscan': 1e6,  # Length of the simulation in unit of round trip
       'f_pmp': 191e12, # Pump Frequency
       'δω_init': 2e9*2*np.pi, # Initial detuning of the pump in rad/s
       'δω_end': -8e9*2*np.pi,  # End detunin of the pump in rad/s
       'μ_sim': [-74,170],  # azimuthal mode to simulate on the left and right side of the pump
       'μ_fit': [-71, 180], # azimuthal mode to fit the dispersion on the left and right side of the pump
        }
```

---

Let’s initialize the class

```python
solver = pyLLE.LLEsovler(sim=sim,
                       res=res,
                       debug=True)
```

## Dispersion Analyse

To plot and retrieve all the data of the dispersion, the method _self.Analyze_ has to be called, resulting in the plot of the integrated dispersion Dint

```python
solver.Analyze(plot=True,
               plottype='all')
```

<iframe frameborder="0" scrolling="no" width="100%" height='400px'src="//plot.ly/~gmoille/34.embed"></iframe>

One can clearly see that because we simulate a larger window with the LLE (orange curve - _LLE simulation_) than the raw data (blue dots), we extrapolate outside of this region. One has to be carefull about this feature as ripple in the integrated dispersion can happen causing zero-crossings which are artefacts

---

A new attribute _disp_ has been created which consists of a dictionary of the different value of the retrieve dispersion

## Temporal Sovler

One can solver the full temporal Lugiato Lefever equation :

$\(t_R \frac{\partial E(t, \tau)}{\partial t} = - \left(\frac{\alpha'}{2} - i\delta_0 \right)E + i \cdot \mathrm{FT}^{-1}\left[ -t_R D_{int}(\omega) \cdot \mathrm{FT}\left[E(t, \tau)\right]\right] + \gamma|E|^2 E + \sqrt{\theta}E_{in}\)$

Test 

<iframe frameborder="0" scrolling="no"  width="100%" height='400px' src="//plot.ly/~gmoille/30.embed"></iframe>