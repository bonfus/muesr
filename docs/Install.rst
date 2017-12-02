Installing Muesr
==================
This section provides an overview and guidance for installing Muesr on
various target platforms.

Prerequisites
-------------
Muesr relies on some 3rd party packages to be fully usable and to
provide you full access to all of its features.

You must have at least the following Python packages installed:

* Python 2.7, 3.1+      (http://www.python.org)
* Numpy 1.6.0+          (http://www.numpy.org)
* mulfc                 (http://github.com/bonfus/muLFC)

Other Python versions or Python implementations might work, but are
(currently) not officially tested or supported.


Additionally, you will need the following packages and libraries to be
installed to use *all* features of Muesr:

========= ========= =============================================== =========================================
Package   Version   Required for                                    Package URL
========= ========= =============================================== =========================================
YAML      >= 2.0.0  :mod:`muesr.i_o.sampleIO`                       http://pyyaml.org/
Spglib    >= 1.6    :mod:`muesr.utilities.symsearch`                http://atztogo.github.io/spglib/
Sympy     >= 1.0    :mod:`muesr.core.magmodel.SMM`                  http://sympy.org
appdirs   >= 1.1    :mod:`muesr.settings`               
XCrysDen  >= 1.0    :mod:`muesr.i_o.xsf.xsf.show_cell` ,            http://www.xcrysden.org
                    :mod:`muesr.i_o.xsf.xsf.show_supercell`          
========= ========= =============================================== =========================================

.. note::
   The muesr distribution ships with a internal version of appdirs which,
   however, may not be up to date.


Compilation from source and system-wide installation
----------------------------------------------------

This is the hard way, but allows to custmize some parts of the installation.
In order to compile the python extension you also need the build tools appropriate
for your system (gcc on Linux, XCode on OS X, Visual Studio or gcc on Windows).
To install the packages you'll need to be superuser.

Use git to clone Muesr and muLFC projects.
First install `muLFC`  ::

    git clone https://github.com/bonfus/muLFC.git
    cd muLFC
    python setup.py install


Next install `muesr`  ::

    git clone https://github.com/bonfus/muesr.git
    cd muesr


In order to install all optional requirements, from within the muesr directory, type :: 

   pip install -r requirements.txt
   
and finally run ::

   python setup.py install


System-wide installation
-------------------------

It is possible to avoid the compilation of the C extension contained
in `lfclib` by directly installing the pre-compiled python wheels.

The wheel can be obtained automatically with pip ::

    pip install --index-url https://testpypi.python.org/pypi mulfc

Finally install Muesr and optional packages ::

    pip install pyyaml spglib
    pip install https://github.com/bonfus/muesr/archive/master.tar.gz


Installation in virtualenv
--------------------------

Virtualenv offers a simple way of virtualizing the Python environment.
This means that you can have a separate collection of python packages 
for running Muesr (and install Muesr itself) without affecting the Python
installation system-wide.

To install Muesr in a virtualenv, first make sure that the command `virtualenv`
is available on your system. If not, please check online what is the 
recommended way of installing virtualenv for your operation system.

To create the virualenv run in a terminal: ::

   virtualenv muesr-env

and to activate the environment (Linux and OsX) ::

   cd muesr-env
   source bin/activate
   
now you can install mulfc and Muesr in the virtualenv with the same commands
reported above ::

    pip install --index-url https://testpypi.python.org/pypi mulfc
    pip install https://github.com/bonfus/muesr/archive/master.tar.gz
   
The `spglib` and `PyYAML` are also recommended ::

   pip install pyyaml spglib


A few notes for Windows users
-----------------------------

In order to install `muesr` on Windows you need a working python environment.
The best user experience is probably provided by Anaconda, which is a
complete Python distribution for scientific data analysis. The following steps assume 
that a working version of `Anaconda <https://www.anaconda.com/download/>`_ is available
on the target system.

Start Anaconda navigator and open an interactive python terminal:

.. image:: anaconda-navigator.png

From within the interactive terminal do: ::

    import pip
    
    pip.main("install --index-url https://testpypi.python.org/pypi mulfc".split())
    pip.main("install spglib".split())
    pip.main("install https://github.com/bonfus/muesr/archive/master.zip".split())



Now you are ready to go! Why not start with a look at the first paragraph
of the :ref:`tutorial` and then move directly to the Muesr :ref:`examples`?




