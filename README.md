# pyLLE ![NIST logo](images/NISTlogo32x32.jpg)

![](https://readthedocs.org/projects/pylle/badge/?version=latest) 
[![](https://img.shields.io/github/license/mashape/apistatus.svg)](licence.txt)

pyLLE is a tool to solver the Lugiato Lefever Equations (LLE)<sup>[1](#ref1)</sup><sup>,</sup><sup>[2](#ref2)</sup><sup>,</sup><sup>[3](#ref3)</sup>in a fast an easy way. Thanks to a user-friendly front-end (and a future UI) in python and a efficient back end in Julia, solving this problem becomes easy and fast. 

For a complete documentation of the package, please visit the [readthedocs page](http://pylle.readthedocs.io/en/latest/index.html)

## Instalation

As pyLLE relies on a Julia back-end, please prior to install this package be sure that Julia is installed on your machine or visit the julia [package downloader page](https://julialang.org/downloads/) to install it. Once Julia installed, the different packages needed to run pyLLE, either python or julia related, will be automatically dowlnloaded and installed 

For a automatic install, just 

```bash
pip install pyLLE
```

For a manual install, download the .zip of the reposotory or clone it and install with the setup.py script 

```bash
git clone https://github.com/gregmoille/pyLLE.git
cd pyLLE
python setup.py install
```


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
       'dispfile': 'h770RW1560.mat' # name of the dispersion file
        }
```



## References

<a name="ref1">1</a>: Luigi A. Lugiato and René Lefever. "Spatial dissipative structures in passive optical systems." Physical review letters 58, no. 21 (1987): 2209.

<a name="ref1">2</a>: Yanne K. Chembo and Curtis R. Menyuk. "Spatiotemporal Lugiato-Lefever formalism for Kerr-comb generation in whispering-gallery-mode resonators." Physical Review A 87, no. 5 (2013): 053852.

<a name="ref1">3</a>: Stéphane Coen, Hamish G. Randle, Thibaut Sylvestre, and Miro Erkintalo. "Modeling of octave-spanning Kerr frequency combs using a generalized mean-field Lugiato–Lefever model." Optics letters 38, no. 1 (2013): 37-39.

<a name="ref1">4</a>: Qing Li, Travis C. Briles, Daron A. Westly, Tara E. Drake, Jordan R. Stone, B. Robert Ilic, Scott A. Diddams, Scott B. Papp, and Kartik Srinivasan. "Stably accessing octave-spanning microresonator frequency combs in the soliton regime." Optica 4, no. 2 (2017): 193-203.