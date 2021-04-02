=========================
Installation Instructions
=========================

pyrtlib can be installed on any computer supporting Python 3.6 (or higher).
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


Create a new virtual environment and activate it:

.. code-block:: console
      
   $ python3 -m venv venv
   $ source venv/bin/activate

Installing from source
----------------------

.. code-block:: console

    $ tar zxvf pyrtlib.tar.gz
    $ cd pyrtlib
    $ python setup.py install --user

pyrtlib is now ready for use from that virtual environment.

.. note::

    To get an up-to-date
    version of pyrtlib, download it directly from `GitHub <https://github.com/slarosa/pyrtlib>`_.