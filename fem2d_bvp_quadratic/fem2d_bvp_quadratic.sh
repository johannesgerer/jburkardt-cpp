#!/bin/bash
#
cp fem2d_bvp_quadratic.hpp /$HOME/include
#
g++ -c -I /$HOME/include fem2d_bvp_quadratic.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling fem2d_bvp_quadratic.cpp"
  exit
fi
#
mv fem2d_bvp_quadratic.o ~/libcpp/$ARCH/fem2d_bvp_quadratic.o
#
echo "Library installed as ~/libcpp/$ARCH/fem2d_bvp_quadratic.o"
