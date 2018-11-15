---
layout: page
title: Example
sidebar_link: true
sidebar_sort_order: 2
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

We now define the simulation parameters. Here we precise a linear detuning ramp of the pump from *δω_init* to *δω_end* relative to the pump mode angular frequency, mode closest to the defined pump frequency *f_pmp*. The simulation length *Tscan* is in unit of round trip, as it is more convenient in the Lugiato-Lefever formalism. It is important to notice that two parameters for the simulation frequency bandwidth have to be defined, *μ_fit* which determined the fit window of the raw data found in *dispfile*, and *μ_sim* which is the number of mode simulated in the LLE, hence could be larger than the fit mode through extrapolation

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

To have a better idea of what happen for a given detuning we can plot the spectra and the temporal profile of the electric field in the resonator



<iframe frameborder="0" scrolling="no" width="100%" height='400px' src="//plot.ly/~gmoille/34.embed"></iframe>

One can clearly see that because we simulate a larger window with the LLE (orange curve - _LLE simulation_) than the raw data (blue dots), we extrapolate outside of this region. One has to be carefull about this feature as ripple in the integrated dispersion can happen causing zero-crossings which are artefacts

---

A new attribute _disp_ has been created which consists of a dictionary of the different value of the retrieve dispersion

## Temporal Sovler

One can solver the full temporal Lugiato Lefever equation 

$$
\begin{align*}
t_R \frac{\partial E(t, \tau)}{\partial t} = &- \left(\frac{\alpha'}{2} - i\delta_0 \right)E +\\ &i \cdot \mathrm{FT}^{-1}\left[ -t_R D_{int}(\omega) \cdot \mathrm{FT}\left[E(t, \tau)\right]\right] +\\& \gamma|E|^2 E + \sqrt{\theta}E_{in}
\end{align*}
$$


The solver implemented in this package is based on a Julia core called through python. Hence to interface both language in an easy way, a hdf5 file is created in a temporary location with all the necessary data to solve the above LLE.

<br>

Hence, to solver the LLE here, we first need to setup the file:

```python
solver.Setup()
```

    -- Solving standard LLE --
            Simulation Parameters
                    R = 23.00 µm
                    Qi = 1.00 M
                    Qc = 1.00 M
                    γ = 1.55
            Simulation Parameters
                    Pin = 150.00 mW
                    Tscan = 1.00 x1e6 Round Trip
                    f_pmp = 191.00 THz
                    δω_init = 2.00 x2π GHz
                    δω_end = -8.00 x2π GHz
                    μ_sim = [-74.00,170.00]
                    μ_fit = [-71.00,180.00]

---

And now we can solve the equation

```python
solver.SolveTemporal()
```

If the user stop the solver (most of the time through a ctrl-c), the Julia solver will be killed hence won't use to much CPU resources

---

When done, we need to retrieve the data

```python
solver.PlotCombPower()
```

## Post Process of Temporal Solver

### Displaying Results

To have a better idea of what happen for a given detuning we can plot the spectra and the temporal profile of the electric field in the resonator

```python
ind = 692
_ = solver.PlotCombSpectra(ind)
```

<!-- PLot All -->
<iframe frameborder="0" scrolling="no" width="100%" height='700px' src="//plot.ly/~gmoille/43.embed"></iframe>

<!-- PLot Sepctra Temporal -->
<iframe frameborder="0" scrolling="no" width="100%" height='400px' src="//plot.ly/~gmoille/46.embed"></iframe>

The temporal profile can also be retrieve, interesting in a case of we are on a soliton state for the given detuning

<!-- PLot Time Temporal -->
<iframe frameborder="0" scrolling="no" width="100%" height='400px' src="//plot.ly/~gmoille/48.embed"></iframe>

### Saving figures

Helper to save the figures, in any format supported by matplotlib have been created. The helper (pylle._llesolver_Latexify) has been designed to make “publication-ready-like” figure, while one is using jupyter notebook/jupyter lab (hence relying on plotly for figure display) or through a python console (hence through matplotlib display). Depending on which type of plot has already been display (comb summary, comb spectra, time profile), this figures will be saved in the specified format

```python
solver.SavePlots2File(basename = 'images/', format = 'pdf')
solver.SavePlots2File(basename = 'images/', format = 'png')
```

### Pickling the Solver

The whole package has been design with simplicity in mind, which mean simplicity to compute but also simplicity to retrieve results. Hence pyLLE is coded in a way that, using the pickle package, one can easily save the state of the solver and retrieve it later, with all the feature presented above still available. It is important to use the self.SaveResults method to pickle the solver as it takes care of dumping the different thread that cannot be pickle in the class.

```python
solver.SaveResults('PickleSolver', path = './')
```

One can reload the state of the simulation:

```python
import pickle as pkl
oldSolver = pkl.load(open('PickleSolver.pkl', 'br'))
```

The different methods can be called with the loaded pickled solver

```python
figly = oldSolver.PlotCombSpectra(692)
```

<!-- PLot Sepctra Temporal -->
<iframe frameborder="0" scrolling="no" width="100%" height='400px' src="//plot.ly/~gmoille/46.embed"></iframe>

### Steady State Solver

One can solver the steady state Lugiato Lefever equation :

$$
\begin{align*}
&- \left(\frac{\alpha'}{2} - i\delta_0 \right)E + i \cdot \mathrm{FT}^{-1}\left[ -t_R D_{int}(\omega) \cdot \mathrm{FT}\left[E(t, \tau)\right]\right] +\\
& \gamma|E|^2 E + \sqrt{\theta}E_{in} = 0
\end{align*}
$$

Using a Newton method to find the root of the equation will result in solving this particular state of the LLE.

The method self.SolveSteadySteate is design for this.

Although, it gives fast results, the accuracy of such solver remains questionable compare to a full temporal resolution of the equation

---

First we need to change the simulation parameters to introduce a fix detuning _δω_

```python
solver.sim['δω'] = -5e9*2*np.pi # more or less what it is as the end of the soliton step
```

```python
steady_fig = solver.SolveSteadySteate()
```
<!-- PLot epctra Steady -->
<iframe frameborder="0" scrolling="no" width="100%" height='400px' src="//plot.ly/~gmoille/50.embed"></iframe>
