#!/bin/bash
#
cp rbf_interp_2d.hpp /$HOME/include
#
g++ -c -g -I/$HOME/include rbf_interp_2d.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling rbf_interp_2d.cpp"
  exit
fi
rm compiler.txt
#
mv rbf_interp_2d.o ~/libcpp/$ARCH/rbf_interp_2d.o
#
echo "Library installed as ~/libcpp/$ARCH/rbf_interp_2d.o"
