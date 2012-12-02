#!/bin/bash
#
cp test_interp_nd.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include test_interp_nd.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_interp_nd.cpp"
  exit
fi
rm compiler.txt
#
mv test_interp_nd.o ~/libcpp/$ARCH/test_interp_nd.o
#
echo "Library installed as ~/libcpp/$ARCH/test_interp_nd.o"
