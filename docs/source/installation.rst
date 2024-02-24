Installation instructions
=========================

pyrtlib can be installed on any computer supporting Python 3.8 (or higher).
The actual installation procedure depends on the operating system. The
instructions below are for Ubuntu and MacOS.

Python Installation (ubuntu)
----------------------------

.. code-block:: console
		
   $ sudo apt update && sudo apt upgrade
   $ sudo apt install python3 python3-pip


Python Installation (macos)
----------------------------

Python3 installation via Homebrew

.. code-block:: console

   $ brew install python3

Python3 can also be installed by downloading the installer from `Python Releases for Mac OS X <https://www.python.org/downloads/mac-osx/>`_.

Installing PyRTlib via PyPi
----------------------------
pyrtlib can be installed via pip from PyPI. To install the package using the following command:

.. code-block:: console
   
   $ pip install pyrtlib

.. note::

    To get an up-to-date
    version of pyrtlib, download it directly from `GitHub <https://github.com/SatCloP/pyrtlib>`_.


Virtual Environment
-------------------

To install virtualenv via pip run:

.. code-block:: console

   $ pip3 install virtualenv


Create a new virtual environment and activate it:

.. code-block:: console

   $ virtualenv -p python3 <desired-path>


Activate the virtualenv:

.. code-block:: console

   $ source <desired-path>/bin/activate


Deactivate the virtualenv:

.. code-block:: console
   
   $ deactivate


Installing PyRTlib from source
------------------------------

Download latest release of pyrtlib source from this `link <https://github.com/SatCloP/pyrtlib/releases/latest>`_.

.. code-block:: console

    $ tar zxvf pyrtlib.tar.gz
    $ cd pyrtlib
    $ <desired-path>/bin/python3 setup.py install

pyrtlib is now ready to be used from that virtual environment. For a quickly test run the following command into the terminal app

.. code-block:: console

   $ <desired-path>/python3 <desired-path>/pyrtlib/hello_spectrum.py

if pyrtlib has been properly installed you should see something like

.. code-block:: console

   $ <desired-path>/python3 <desired-path>/pyrtlib/hello_spectrum.py
   Progress: |██████████████████████████████████████████████████| 100.0% Complete
   Hello Spectrum!

               tbtotal  tbatm         tmr  tmrcld     tauwet    taudry  tauliq  tauice
   18.7000   298.689123    0.0  286.716080     0.0   0.069040  0.012013     0.0     0.0
   23.8000   297.014923    0.0  286.634107     0.0   0.214403  0.015643     0.0     0.0
   31.4000   298.285354    0.0  285.140186     0.0   0.076330  0.025881     0.0     0.0
   50.3000   290.594440    0.0  274.191598     0.0   0.124585  0.316968     0.0     0.0
   52.6100   278.442378    0.0  267.163248     0.0   0.134824  0.924593     0.0     0.0
   53.2400   270.032638    0.0  262.487813     0.0   0.137720  1.458056     0.0     0.0
   53.7500   259.296109    0.0  255.080703     0.0   0.140096  2.219325     0.0     0.0
   89.0000   295.336793    0.0  286.913337     0.0   0.370017  0.047366     0.0     0.0
   115.5503  283.409636    0.0  274.910320     0.0   0.634700  0.435743     0.0     0.0
   116.6503  273.105313    0.0  265.583070     0.0   0.647756  0.864176     0.0     0.0
   117.3503  258.382394    0.0  253.279983     0.0   0.656168  1.551855     0.0     0.0
   117.5503  251.887074    0.0  247.840191     0.0   0.658587  1.892017     0.0     0.0
   119.9503  252.319901    0.0  248.289379     0.0   0.688148  1.857808     0.0     0.0
   120.1503  258.829337    0.0  253.792452     0.0   0.690658  1.519190     0.0     0.0
   120.8503  273.470564    0.0  266.281272     0.0   0.699499  0.837028     0.0     0.0
   121.9503  283.508571    0.0  275.765375     0.0   0.713579  0.414934     0.0     0.0
   164.7750  287.382258    0.0  285.293882     0.0   1.912160  0.019109     0.0     0.0
   166.2250  286.768856    0.0  284.923583     0.0   2.061262  0.019146     0.0     0.0
   174.9100  279.316272    0.0  279.136791     0.0   4.721552  0.019642     0.0     0.0
   177.2100  274.918510    0.0  274.902966     0.0   7.354952  0.019836     0.0     0.0
   178.4100  271.637064    0.0  271.635743     0.0   9.944304  0.019946     0.0     0.0
   179.9100  265.916650    0.0  265.916645     0.0  15.761551  0.020091     0.0     0.0
   181.3100  258.183942    0.0  258.183942     0.0  26.052880  0.020233     0.0     0.0
   185.3100  258.248076    0.0  258.248076     0.0  26.149293  0.020672     0.0     0.0
   186.7100  265.558982    0.0  265.558979     0.0  16.344414  0.020837     0.0     0.0
   188.2100  270.889844    0.0  270.889228     0.0  10.732092  0.021020     0.0     0.0
   189.4100  273.904425    0.0  273.897462     0.0   8.196430  0.021170     0.0     0.0
   191.7100  277.820891    0.0  277.740367     0.0   5.586945  0.021468     0.0     0.0

   PyRTlib successfully installed


Build and run the Docker image
===============================

To build docker image it is necessary to download the latest  pyrtlib release from this `link <https://github.com/SatCloP/pyrtlib/releases/latest>`_ and then run the following command from you prefer terminal.

.. code-block:: console

   $ tar zxvf pyrtlib.tar.gz
   $ cd pyrtlib

From within pyrtlib folder run the following docker command to build the docker image

.. code-block:: console

   $ docker build --pull --rm -f "Dockerfile" -t pyrtlib:latest "." 
   $ docker run --rm -it  pyrtlib:latest

To test run the example script from within the docker image

.. code-block:: console

   $ root@993587e5fea9:/home/dev/pyrtlib# python3 hello_spectrum.py

My first run with PyRTlib (Colab Notebook)
==========================================

To run the example script in a Google Colab Notebook, you can use the following code:

.. code-block:: console

   !pip install pyrtlib
   !python3 hello_spectrum.py

.. note::

    The example script is available at `this link <https://colab.research.google.com/github/SatCloP/pyrtlib/blob/main/docs/source/notebook/first_run.ipynb>`_.

.. code-block:: console

   !wget https://raw.githubusercontent.com/SatCloP/pyrtlib/main/pyrtlib/hello_spectrum.py
   !python3 hello_spectrum.py

