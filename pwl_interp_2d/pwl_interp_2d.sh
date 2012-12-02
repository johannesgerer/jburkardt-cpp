#!/bin/bash
#
cp pwl_interp_2d.hpp /$HOME/include
#
g++ -c -g -I/$HOME/include pwl_interp_2d.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pwl_interp_2d.cpp"
  exit
fi
rm compiler.txt
#
mv pwl_interp_2d.o ~/libcpp/$ARCH/pwl_interp_2d.o
#
echo "Library installed as ~/libcpp/$ARCH/pwl_interp_2d.o"
