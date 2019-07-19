---
layout: page
title: Installation
sidebar_link: true
sidebar_sort_order: 1
---

## Requirement

As pyLLE relies on a Julia back-end, please prior to installing this package be sure that Julia is installed on your machine or visit the julia [package download page](https://julialang.org/downloads/) to install it by selecting &#9888; **v1.1.1** &#9888;.

**Windows users**: Please, keep julia in the default directory during the installation (i.e. ~\AppData\Local\Julia-1.1.1\ for windows).

**Mac Os User**: You would need to add the julia binary to the path. The easiest way to do it is to create a simlink in the terminal


## Auto installation


After installing Julia on the machine just run:

```
pip install pyLLE
```

## Manual Installation


Download the .zip [repository](https://github.com/gregmoille/pyLLE/archive/master.zip) or clone it using

```
git clone https://github.com/gregmoille/pyLLE.git
cd pyLLE
python setup.py install
```
Just a heads up, the installation can vary in time, especially because of Julia that might rebuild the cache. If the Julia location is custom, please before running the last line, change in the *setup.py* line 18 to the correct location, as in *pyLLE/_llesolver.py* line **508** to point to the correct location. Thanks

## Python Dependencies

- scipy
- numpy
- matplotlib
- h5py
- prettytable
- matplotlib
- ipdb

## Julia Dependencies

- HDF5
- FFTW
- Base
- LinearAlgebra

## Checking the Installation

### Checking Python Installation

- Open a python command prompt. You can simply type python in a terminal (Linux/Mac Os) or a cmd window (Windows), or you can launch spyder, jupyter or a ipython instance.
- Inside python, try to import the different package that pyLLE requires to check if they were installed correctly:

```python
import scipy
import plotly
import numpy
import matplotlib
import h5py
import prettytable
import matplotlib
import ipdb
```

If any errors show up due to a missing package please:

- Open a terminal (Linux/MacOs) or cmd window (Windows). **The following commands should not be type in a python instance** but directly in the shell.
- Install the missing package by typing:

```
pip install <missing package name>
```

### Checking Julia Installation

- open a Julia shell
- type the following commands:

```julia
using Pkg

Pkg.add("HDF5")
Pkg.update("HDF5")

Pkg.add("FFTW")
Pkg.update("FFTW")

Pkg.add("LinearAlgebra")
Pkg.update("LinearAlgebra")
```

- You should not have any errors at that level usually.
- If you have an issue with the installation of a Julia package, directly enquire to the faulty pacakge through their github page as they might have solution to install it (the odd that is happens is really small though)
