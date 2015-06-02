#!/bin/bash
#
cp hypercube_integrals.hpp /$HOME/include
#
g++ -c -I /$HOME/include hypercube_integrals.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling hypercube_integrals.cpp"
  exit
fi
#
mv hypercube_integrals.o ~/libcpp/$ARCH/hypercube_integrals.o
#
echo "Library installed as ~/libcpp/$ARCH/hypercube_integrals.o"
