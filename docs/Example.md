---
layout: page
title: Example
sidebar_link: true
sidebar_sort_order: 2
---

This example is accessible through a Jupyter notebook available in the [example folder](https://github.com/gregmoille/pyLLE/tree/master/example). The use of Jupyter based on the convenience of its interactivity, but obviously one can script pyLLE and use an ipython console or a python bash command, where one would use *matplotlib* instead of *plotly*. Check the [documentation of the module](https://gregmoille.github.io/pyLLE/pyLLE.html) for the difference in output in both cases.

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

We now define the simulation parameters. Here we specify a linear detuning ramp of the pump from *δω_init* to *δω_end* relative to the pump mode angular frequency, which is the mode closest to the defined pump frequency *f_pmp*. The simulation length *Tscan* is in units of round trip, as it is more convenient in the Lugiato-Lefever formalism. It is important to notice that two parameters for the simulation frequency bandwidth have to be defined, *μ_fit* which determines the fit window of the raw data found in *dispfile*, and *μ_sim* which is the number of modes simulated in the LLE. *μ_fit* can be larger than *μ_sim* through extrapolation.

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
solver = pyLLE.LLEsolver(sim=sim,
                       res=res,
                       debug=True)
```

## Dispersion Analysis

To plot and retrieve all the dispersion data, the method _self.Analyze_ has to be called, resulting in the plot of the integrated dispersion Dint

```python
solver.Analyze(plot=True,
               plottype='all')
```

<iframe frameborder="0" scrolling="no" width="100%" height='400px' src="//plot.ly/~gmoille/34.embed"></iframe>

Here, we extrapolate (orange curve - _LLE simulation_) outside of the spectral window containing the raw dispersion data (blue dots), to enable a simulation over a broader spectral range. One has to be carefull about using this extrapolation feature, as ripples in the integrated dispersion can result, creating zero-crossings which are artefacts that will influence the LLE simulation results. 

---

A new attribute _disp_ has been created, and consists of a dictionary of the different values of the retrieve dispersion.

## Temporal Sovler

One can solver the full temporal Lugiato Lefever equation 

$$
\begin{align*}
t_R \frac{\partial E(t, \tau)}{\partial t} = &- \left(\frac{\alpha'}{2} - i\delta_0 \right)E + \\ &i \cdot \mathrm{FT}^{-1}\left[ -t_R D_{int}(\omega) \cdot \mathrm{FT}\left[E(t, \tau)\right]\right] +\\& \gamma|E|^2 E + \sqrt{\theta}E_{in}
\end{align*}
$$


The solver implemented in this package is based on a Julia core called through python. Hence, to interface to both languages in an easy way, a hdf5 file is created in a temporary location with all the necessary data to solve the above LLE.

<br>

To solve the LLE, we first need to set up the file:

```python
solver.Setup()
```

And now we can solve the equation

```python
solver.SolveTemporal()
```

If the user stops the solver (most of the time through a ctrl-c command), the Julia solver will be killed, and hence won't use too much CPU resources.

---

When done, we need to retrieve the data:

```python
solver.PlotCombPower()
```

Here, we have an overall indication of the behavior of the system through a 2D plot of the comb power vs. frequency and LLE step (top), a 2D plot of the comb power vs. fast time and LLE step (middle), and a 1D plot of the comb power vs. LLE step (bottom).  

<!-- PLot All -->
<iframe frameborder="0" scrolling="no" width="100%" height='700px' src="//plot.ly/~gmoille/43.embed"></iframe>


## Post Processing of the Temporal Solver

### Displaying Results

To have a better idea of what is happening for a given laser-cavity mode detuning, we can plot the spectrum and the temporal profile of the electric field in the resonator at a fixed LLE step:

```python
ind = 570
_ = solver.PlotCombSpectra(ind)
```

<!-- PLot Sepctra Temporal -->
<iframe frameborder="0" scrolling="no" width="100%" height='400px' src="//plot.ly/~gmoille/46.embed"></iframe>

The temporal profile can also be retrieved, which can be interesting, for example, in the case when we are on a soliton state for the given detuning:

<!-- PLot Time Temporal -->
<iframe frameborder="0" scrolling="no" width="100%" height='400px' src="//plot.ly/~gmoille/48.embed"></iframe>

### Saving figures

Helper routines to save the figures in any format supported by matplotlib have been created. The helper (pylle._llesolver_Latexify) has been designed to make “publication-ready-like” figures, while one is using jupyter notebook/jupyter lab (hence relying on plotly for figure display) or through a python console (hence through matplotlib display). Depending on which type of plot has already been displayed (comb summary, comb spectra, time profile), these figures will be saved in the specified format.

```python
solver.SavePlots2File(basename = 'images/', format = 'pdf')
solver.SavePlots2File(basename = 'images/', format = 'png')
```

### Pickling the Solver

The whole package has been design with simplicity in mind, which means not only simplicity to compute, but also simplicity to retrieve results. Hence pyLLE is coded in a way that, using the pickle package, one can easily save the state of the solver and retrieve it later, with all the features presented above still available. It is important to use the self.SaveResults method to pickle the solver as it takes care of dumping the different threads that cannot be pickled in the class.

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
figly = oldSolver.PlotCombSpectra(570)
```

<!-- PLot Sepctra Temporal -->
<iframe frameborder="0" scrolling="no" width="100%" height='400px' src="//plot.ly/~gmoille/46.embed"></iframe>

### Steady State Solver

One can solve the steady state Lugiato Lefever equation :

$$
\begin{align*}
&- \left(\frac{\alpha'}{2} - i\delta_0 \right)E + i \cdot \mathrm{FT}^{-1}\left[ -t_R D_{int}(\omega) \cdot \mathrm{FT}\left[E(t, \tau)\right]\right] +\\
& \gamma|E|^2 E + \sqrt{\theta}E_{in} = 0
\end{align*}
$$

Using a Newton's method to find the root of the equation will result in solving this particular state of the LLE.

The method self.SolveSteadySteate is designed for this. Although, it gives fast results, the accuracy of such a solver remains questionable compared to a full temporal resolution of the equation (in particular, the convergence of the Newton's method is not always assured). 

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
