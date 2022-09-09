# pyLLE ![NIST logo](images/NISTlogo32x32.jpg)
[![](https://img.shields.io/static/v1.svg?label=docs&message=passing&color=green&style=flat)](https://gregmoille.github.io/pyLLE/)
![](https://img.shields.io/static/v1.svg?label=version&message=4.0.0&color=9cf&style=flat)
[![](https://img.shields.io/static/v1.svg?label=DOI&message=10.6028/jres.124.012&color=blue&style=flat)](https://doi.org/10.6028/jres.124.012)

## What's new

Quite of a new version, v4.0.0 introduce a lot of new stuff including: 
- Complete recoding of the core with better codding and commenting which hopefully makes it easier to implement custom features 
- Modification of interfacing with parameters and results through method attributes instead of bulky dictionaries
- Stability of the half step Fourier method, allowing to use a soliton solution as an original state for the LLE 
- Allowing arbitrary number of driving pump, according to Taheri et al. and our paper on Nature Communication (Moille et al. )
- Julia compatibility with version 1.1 and above 

## How to Cite Us?

Please, if you use this package and it helps you with your research adn publication, cite us in your paper.Not only it allows us to have a better idea of new things people are interested in and how to keep improving the solver, but it also help us building a community where every body could help maintaining the solver to suit better the needs of everybody.
You can cite our paper published in the Journal of Research of National Institute of Standards and Technology available [here](https://doi.org/10.6028/jres.124.012), with the following bibtex entry:

```latex
@article{moille_pyLLE,
      author = {Gregory Moille and Qing Li and Xiyuan Lu and Kartik Srinivasan},
      title = {pyLLE: a Fast and User Friendly Lugiato-Lefever Equation Solver},
      year = {2019},
      volume = {124},
      pages = {124012},
      month = {2019-05-24},
      journal = {Journal of Research of NIST},
       doi = {https://doi.org/10.6028/jres.124.012},
     }
```

## Example [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/gregmoille/pyLLE/HEAD?labpath=example%2FTemporalDualPump.ipynb)

As pyLLE relies on a Julia back-end, please prior to installing this package be sure that Julia is installed on your machine or visit the julia [package download page](https://julialang.org/downloads/) to install it. The code should now work with any recent version of Julia.


**Windows users**: Please, keep julia in the default directory during the installation (i.e. ~\AppData\Local\Julia-1.1.1\ for windows).

**Mac Os User**: You would need to add the julia binary to the path. The easiest way to do it is to create a simlink in the terminal

```bash
ln -s /Applications/Julia-<version>.app/Contents/Resources/julia/bin/julia /usr/local/bin/julia
```

Once Julia installed, the different packages needed to run pyLLE, either python or julia related, will be automatically downloaded and installed. Just a heads up, the installation of the package can vary in time, especially because of Julia that might rebuild the cache.
For a automatic install, just pip it :

```bash
pip install pyLLE
```
For a manual install, download the .zip of the repository or clone it and install with the setup.py script


If the julia location is custom, please before installing change in the setup.py, line 18 to the correct location, as in pyLLE/llesolver.py line 430 to point to the correct location. Thanks


## Checking that everything works correctly

Launch a julia console and within type the commands:

```julia
using HDF5
using FFTW
using LinearAlgebra
```


if any of the previous command throw an issue, mostly it is because it is not installed. One way to fix it is to remove the installed packaged to remove the cache

- for linux and mac os user: remove everything in ~/.julia/
- for windows users: remove everything in C:\Users\<your user name>\.julia\

Then enter the pacakge manager for julia by typing in the julia console:

```julia
julia>]
```

then
```julia
(v1.1) pkg>add HDF5
(v1.1) pkg>add FFTW
```

## Example [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/gregmoille/pyLLE/HEAD?labpath=example%2FTemporalDualPump.ipynb)

A complete example is available in the example directory [notebook](https://github.com/gregmoille/pyLLE/blob/master/example/TemporalDualPump.ipynb) with the corresponding file needed in the folder. 
You can also access the [interactive binder](https://mybinder.org/v2/gh/gregmoille/pyLLE/HEAD?labpath=example%2FTemporalDualPump.ipynb) example to have a better idea of what's going on: 

## Works Using pyLLE 

If you want to be featured here, shoot me an email! I try to keep it up to date but this is not a priority, yet I would love to hear anybody who uses it!

- Gregory Moille, Xiyuan Lu, Ashutosh Rao, Qing Li, Daron A. Westly, Leonardo Ranzani, Scott B. Papp, Mohammad Soltani, and Kartik Srinivasan "_Kerr-Microresonator Soliton Frequency Combs at Cryogenic Temperatures_," Phys. Rev. Applied 12, 034057 (2019)
- Gregory Moille, Qing Li, Travis C. Briles, Su-Peng Yu, Tara Drake, Xiyuan Lu, Ashutosh Rao, Daron Westly, Scott B. Papp, and Kartik Srinivasan "_Broadband resonator-waveguide coupling for efficient extraction of octave-spanning microcombs_," Optics Letters Vol. 44, Issue 19, pp. 4737-4740 (2019)
- Lin Chang, Weiqiang Xie, Haowen Shu, Qifan Yang, Boqiang Shen, Andreas Boes, Jon D. Peters, Warren Jin, Songtao Liu, Gregory Moille, Su-Peng Yu, Xingjun Wang, Kartik Srinivasan, Scott B. Papp, Kerry Vahala, John E. Bowers "_Ultra-efficient frequency comb generation in AlGaAs-on-insulator microresonators_," arXiv:1909.09778 (2019)
- Schuttrups, B. (2020). "_Modelling nonlinear optical pulse propagation using pseudo-spectral methods_" (Master's thesis, University of Twente).
- Moille, G., Chang, L., Xie, W., Rao, A., Lu, X., Davanco, M. _et al._  "_Dissipative Kerr Solitons in a III‐V Microresonator_". Laser & Photonics Reviews, 14(8), 2000022 (2020)
- Weng, H., Liu, J., Afridi, A. A., Li, J., Dai, J., Ma, X. _et al._ "_Octave-spanning Kerr frequency comb generation with stimulated Raman scattering in an AlN microresonator_". Optics Letters, 46(3), 540-543. (2021)
- Weng, H., Liu, J., Afridi, A. A., Li, J., Dai, J., Ma, X. _et al._ "_Directly accessing octave-spanning dissipative Kerr soliton frequency combs in an AlN microring resonator_" Photonics Research (2021)
- Moille, G., Westly, D., Orji, N. G., & Srinivasan, K. "_Tailoring broadband Kerr soliton microcombs via post-fabrication tuning of the geometric dispersion_". Applied Physics Letters, 119(12), 121103 (2021)
- Weng, H., Liu, J., Afridi, A. A., Li, J., Dai, J., Zhang, Y., _et al._ "_Perfect soliton crystal in an AlN microresonator"_. In CLEO: QELS_Fundamental Science (pp. JTh3A-31). Optical Society of America. (2021)
- Moille, G., Westly, D., Simelgor, G., & Srinivasan, K. "_Impact of the precursor gas ratio on dispersion engineering of broadband silicon nitride microresonator frequency combs_". Optics Letters, 46(23), 5970-5973 (2021).
