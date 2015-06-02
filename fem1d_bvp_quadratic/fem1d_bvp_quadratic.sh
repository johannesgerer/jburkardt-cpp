#!/bin/bash
#
cp fem1d_bvp_quadratic.hpp /$HOME/include
#
g++ -c -I /$HOME/include fem1d_bvp_quadratic.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling fem1d_bvp_quadratic.cpp"
  exit
fi
#
mv fem1d_bvp_quadratic.o ~/libcpp/$ARCH/fem1d_bvp_quadratic.o
#
echo "Library installed as ~/libcpp/$ARCH/fem1d_bvp_quadratic.o"
