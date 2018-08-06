Installation
==================

Requirement
--------------------

pyLLE relies on a **python >3.4** front-end for a user friendly interface – either through script commands or GUI - and on a **julia** back-end for fast computation. Hence, one need to `download <https://julialang.org/downloads/>`_ and install julia first before running the installation of pyLLE. 

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
