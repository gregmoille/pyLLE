Installation
==================

Requirement
--------------------

pyLLE relies on a **python >3.4** front-end for a user friendly interface – either through script commands or GUI - and on a **julia v0.6.4** back-end for fast computation. Hence, one need to `download <https://julialang.org/downloads/>`_ and install julia first before running the installation of pyLLE. Please, keep julia in the default directory during the installation (i.e. ~\AppData\Local\Julia-0.6.4\). If not, please go to the manual installation

The different julia packages required for solving the LLE and interfacing with python will be installed while pyLLE is set up. 


Auto installation
--------------------

After installing julia on the machine just run: 

.. code:: bash

    pip install pyLLE


Manual Installation
--------------------

Download the .zip `repository <https://github.com/gregmoille/pyLLE/archive/master.zip>`_ or clone it using 

.. code:: bash

    git clone https://github.com/gregmoille/pyLLE.git
    cd pyLLE
    python setup.py install

Just a heads up, the installation can vary in time, especially because of Julia that might rebuilds the cache. If the julia location is custom, please before running the last line, change in the setup.py line 18 to the correct location, as in pyLLE/llesolver.py line 430 to point to the correct location. Thanks

Python Dependencies
--------------------

- scipy
- numpy
- matplotlib
- h5py
- prettytable
- matplotlib

Julia Dependencies
--------------------

- HDF5
- ProgressMeter