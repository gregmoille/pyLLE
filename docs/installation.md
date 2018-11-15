---
layout: page
title: Installation
sidebar_link: true
sidebar_sort_order: 1
---

## Requirement


pyLLE relies on a **python >3.4** front-end for a user friendly interface – either through script commands or GUI - and on a **Julia v0.6.4** back-end for fast computation. Hence, one need to [download](https://Julialang.org/downloads/oldreleases.html) and install Julia first before running the installation of pyLLE. **It is important to keep version 0.6.4 as there is some issue with HDF5 with further released**. Please, keep Julia in the default directory during the installation (i.e. ~\AppData\Local\Julia-0.6.4\). If not, please go to the manual installation

The different Julia packages required for solving the LLE and interfacing with python will be installed while pyLLE is set up. 


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
Just a heads up, the installation can vary in time, especially because of Julia that might rebuilds the cache. If the Julia location is custom, please before running the last line, change in the _setup.py_ line 18 to the correct location, as in _pyLLE/_llesolver.py_ line **508** to point to the correct location. Thanks

## Python Dependencies

- scipy
- numpy
- matplotlib
- h5py
- prettytable
- matplotlib

## Julia Dependencies

- HDF5
- ProgressMeter