#!/bin/bash
#
g++ -c -g fem2d_poisson.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem2d_poisson.cpp"
  exit
fi
rm compiler.txt
#
mv fem2d_poisson.o ~/libcpp/$ARCH
#
echo "Partial program installed as ~/libcpp/$ARCH/fem2d_poisson.o"
