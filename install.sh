#!/bin/bash

#compile the symlib submodule
cd enumlib/symlib/src
export F90=gfortran
make

#make the enumeration library itself
cd ../../src
make
make enum.x
make makestr.x
