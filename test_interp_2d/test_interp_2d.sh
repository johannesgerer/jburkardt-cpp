#!/bin/bash
#
cp test_interp_2d.hpp /$HOME/include
#
g++ -c -I /$HOME/include test_interp_2d.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling test_interp_2d.cpp."
  exit
fi
#
mv test_interp_2d.o ~/libcpp/$ARCH/test_interp_2d.o
#
echo "Library installed as ~/libcpp/$ARCH/test_interp_2d.o"
