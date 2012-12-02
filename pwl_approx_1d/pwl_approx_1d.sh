#!/bin/bash
#
cp pwl_approx_1d.hpp /$HOME/include
#
g++ -c -g -I/$HOME/include pwl_approx_1d.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pwl_approx_1d.cpp"
  exit
fi
rm compiler.txt
#
mv pwl_approx_1d.o ~/libcpp/$ARCH/pwl_approx_1d.o
#
echo "Library installed as ~/libcpp/$ARCH/pwl_approx_1d.o"
