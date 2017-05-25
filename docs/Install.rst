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

Other Python versions or Python implementations might work, but are
(currently) not officially tested or supported by the Muesr
distribution.

.. note::
   Windows is likely to be supported with minor changes to the code.

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


Compilation and installation
----------------------------

This is the hard way. Make your git project muesr directory, move into it and use ::

  git clone https://github.com/bonfus/muesr.git

In order to compile the python extension you also need the build tools appropriate
for your system (gcc on Linux, XCode on OS X, Visual Studio or gcc on Windows).

If you have all the prerequisites and you are in the muesr directory, type:: 

   make clean
   make all
   make install

If you do not want to compile, you can use instead the following `wheels: <https://packaging.python.org/wheel_egg/>`_

Direct installation
-------------------
Use the python way of installing the package. First, `download it <https://github.com/bonfus/muesr/archive/master.zip>`_ and unzip. Move to the muesr directory and simply type ::
   
   python setup.py test
   python setup.py install

You can also use pip ::

   pip install muesr-0.1.tar.gz

(you'll probably need to be superuser)
