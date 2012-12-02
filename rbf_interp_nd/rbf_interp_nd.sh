#!/bin/bash
#
cp rbf_interp_nd.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include rbf_interp_nd.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling rbf_interp_nd.cpp"
  exit
fi
rm compiler.txt
#
mv rbf_interp_nd.o ~/libcpp/$ARCH/rbf_interp_nd.o
#
echo "Library installed as ~/libcpp/$ARCH/rbf_interp_nd.o"
