#!/bin/bash
#
cp shepard_interp_nd.hpp /$HOME/include
#
g++ -c -g -I/$HOME/include shepard_interp_nd.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling shepard_interp_nd.cpp"
  exit
fi
rm compiler.txt
#
mv shepard_interp_nd.o ~/libcpp/$ARCH/shepard_interp_nd.o
#
echo "Library installed as ~/libcpp/$ARCH/shepard_interp_nd.o"
