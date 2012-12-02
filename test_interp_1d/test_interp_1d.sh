#!/bin/bash
#
cp test_interp_1d.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include test_interp_1d.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_interp_1d.cpp"
  exit
fi
rm compiler.txt
#
mv test_interp_1d.o ~/libcpp/$ARCH/test_interp_1d.o
#
echo "Library installed as ~/libcpp/$ARCH/test_interp_1d.o"
