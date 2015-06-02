#!/bin/bash
#
cp multigrid_poisson_1d.hpp /$HOME/include
#
g++ -c multigrid_poisson_1d.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling multigrid_poisson_1d.cpp."
  exit
fi
#
mv multigrid_poisson_1d.o ~/libcpp/$ARCH/multigrid_poisson_1d.o
#
echo "Library installed as ~/libcpp/$ARCH/multigrid_poisson_1d.o"
