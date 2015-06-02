#!/bin/bash
#
cp jacobi_eigenvalue.hpp /$HOME/include
#
g++ -c -I/$HOME/include jacobi_eigenvalue.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling jacobi_eigenvalue.cpp"
  exit
fi
#
mv jacobi_eigenvalue.o ~/libcpp/$ARCH/jacobi_eigenvalue.o
#
echo "Library installed as ~/libcpp/$ARCH/jacobi_eigenvalue.o"
