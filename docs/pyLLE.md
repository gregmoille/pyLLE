---
layout: page
title: pyLLE Modules
sidebar_link: true
sidebar_sort_order: 6
---

# _llesolver

## MyLogger
```python
MyLogger(self, fname)
```

Custom made logger as the logger default package cannot be pickled

## Latexify
```python
Latexify(self, **kwarg)
```

Class that handle saving the figures in a nice way compatible with
the column/page size of different latex template

- input [] = optional:
    - figname = name to save the figure (without extension)
    - fig = matplotlib handle to the figure
    - [fig_width]: default = 1column
    - [frmt]: default = pdf
    - [fig_height] : default = 6.5
    - [font_size] : default = 8

## LLEsovler
```python
LLEsovler(self, **kwargs)
```

Class to solve the Lugiato Lefever Equation
Initialization input ([]=facultative):

- res \<dict\>
    - Qi \<float\>: intrinsic Q of the resonator
    - Qc \<float\>: coupling Q of the resonator
    - R \<float\>: ring radius
    - gamma \<float\>: Non linear index of the material
    - dispfile \<str\> : str pointing to a .csv file where the azimuthal mode orders and corresponding resonances are saved
- sim \<dict\>
    - Tscan \<float\>: length of the simulation (in unit of round trip)
    - mu_fit \<list\>: number of mode to fit
    - mu_sim \<list\>: number of mode to simulate
    - domega_init \<float\>: initial detuning of the pump
    - domega_end \<float\>: final detuning of the pump
    - [domga_stop] \<float\>: where to stop the scan in detuning but keep doing the simulation
- debug \<bool\>: Save a trace in a logfile in the working directory of the different actions pyLLE perform (default = True)

# LLEsovler
```python
LLEsovler(self, **kwargs)
```

Class to solve the Lugiato Lefever Equation
Initialization input ([]=facultative):

- res \<dict\>
    - Qi \<float\>: intrinsic Q of the resonator
    - Qc \<float\>: coupling Q of the resonator
    - R \<float\>: ring radius
    - gamma \<float\>: Non linear index of the material
    - dispfile \<str\> : str pointing to a .csv file where the azimuthal mode orders and corresponding resonances are saved

- sim \<dict\>
    - Tscan \<float\>: length of the simulation (in unit of round trip)
    - mu_fit \<list\>: number of mode to fit
    - mu_sim \<list\>: number of mode to simulate
    - domega_init \<float\>: initial detuning of the pump
    - domega_end \<float\>: final detuning of the pump
    - [domga_stop] \<float\>: where to stop the scan in detuning but keep doing the simulation

- debug \<bool\>: Save a trace in a logfile in the working directory of the different actions pyLLE perform (default = True)

### hbar
float(x) -\> floating point number

Convert a string or number to a floating point number, if possible.
### Analyze
```python
LLEsovler.Analyze(self, plot=False, f=None, ax=None, label=None, plottype='all', zero_lines=True, mu_sim=None)
```

Call pyLLE.analyzedisp.AnalyzeDisp to get the dispersion of the resonator we want to simulate

### Setup
```python
LLEsovler.Setup(self)
```

Setup the simulation for the Julia back-end.
Save the two main dictionary self.sim and self.res into a readable hdf5 file for Julia in the temporary location define by the os

### SolveTemporal
```python
LLEsovler.SolveTemporal(self, tol=0.001, maxiter=6, step_factor=0.1)
```

Call Julia to solve the LLE

### SolveSteadySteate
```python
LLEsovler.SolveSteadySteate(self)
```

Newton Method to find the root of the steady state equation

### RetrieveData
```python
LLEsovler.RetrieveData(self)
```

Load the output hdf5 saved by julia and transform it in a user-friendly dictionary to be more pythonistic

### PlotCombPower
```python
LLEsovler.PlotCombPower(self, do_matplotlib=False)
```

Plot a figure with 3 subplots.

- Top subplot = map of the spectra for the steps taken by the LLE (step sub-sampled to be 1000)
- middle subplot = temporal map of the intensity inside the resonator for the steps of the LLE
- bottom subplot = normalized comb power

- Output
    - f, ax:  handle of figure and axes of the matplotlib figure displayed

### PlotCombSpectra
```python
LLEsovler.PlotCombSpectra(self, ind, f=None, ax=None, label=None, pwr='both', do_matplotlib=False, plot=True)
```

Plot the spectra for a given index in the 1000 sub-sampled LLE steps

- Input
    - ind \<ind\>: index in the LLE step to plot the spectra
    - f \<obj\>:  matplotlib figure handle (if None, new figure)
    - ax \<obj\>: matplotlib axe handle
    - label \<str\>: label for the legend
    - pwr \<str\>: 'both', 'ring', 'wg' depending on the spectra wanted (inside the ring, the waveguide or both)
- Output
    - freq \<numpy.array\>: frequency in Hz
    - Sout \<numpy.array\>: spectral density of power in the waveguide (dBm)
    - Sring \<numpy.array\>: spectral density of power in the ring (dBm)
    - f \<obj\>:  matplotlib figure handle
    - ax \<obj\>: matplotlib axes handle

### PlotSolitonTime
```python
LLEsovler.PlotSolitonTime(self, ind, f=None, ax=None, label=None, do_matplotlib=False)
```

Plot the spectra for a given index in the 1000 sub-sampled LLE step

**Input**

    - ind \<ind\>: index in the LLE step to plot the spectra
    - f \<obj\>:  matplotlib figure handle (if None, new figure)
    - ax \<obj\>: matplotlib axe handle
    - label \<str\>: label for the legend

**Output**

    - τ \<obj\>: Time in the resonator
    - U \<numpy.array\>: Temporal Electric field for the given step of the LLE
    - f \<obj\>: matplotlib figure handle
    - ax \<obj\>: matplotlib axe handle

### SaveResults
```python
LLEsovler.SaveResults(self, fname, path='./')
```

Save the whole class with pickle to be able to easilly call it back or retrieve the results after saving

**Input**

    - fname \<str\>: name to save. The '.pkl' extension will be added
    - path \<str\>: path to save the results (defaults './')

