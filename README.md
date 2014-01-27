Molecule parametrization
========================

This repository consists of the following projects.

* Symmetrization
* Charge-group partitioning
* Maximal common fragment identification

Dependencies
------------

* LEMON 1.3
* Qt 5.x (only for cgp-demo)

Compilation instructions
------------------------

First, LEMON 1.3 needs to be installed:

    wget http://lemon.cs.elte.hu/pub/sources/lemon-1.3.tar.gz
    tar xvzf lemon-1.3.tar.gz
    cd lemon-1.3
    cmake -DCMAKE_INSTALL_PREFIX=~/lemon
    make install
    
Note: On Mac OS 10.9, comment out the following two lines in `CMakeLists.txt` before doing `make install`

    #ADD_SUBDIRECTORY(demo) 
    #ADD_SUBDIRECTORY(tools)
    
You can remove the LEMON sources now, i.e., `rm -rf lemon-1.3`. 

To compile the projects:

    mkdir build
    cd build
    cmake ..
    make
    
In case auto-detection of LEMON fails, do

    cmake -DLIBLEMON_ROOT=~/lemon ..
    
To compile cgp-demo:

    cd src/cgp-demo
    qmake
    make
    
Running
-------

To run `fragments`

    ./fragments -no-json -s 1 ../../../data/CGP/fragments/16841.lgf ../../../data/CGP/fragments/16842.lgf
