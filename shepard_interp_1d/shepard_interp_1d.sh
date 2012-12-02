#!/bin/bash
#
cp shepard_interp_1d.hpp /$HOME/include
#
g++ -c -g -I/$HOME/include shepard_interp_1d.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling shepard_interp_1d.cpp"
  exit
fi
rm compiler.txt
#
mv shepard_interp_1d.o ~/libcpp/$ARCH/shepard_interp_1d.o
#
echo "Library installed as ~/libcpp/$ARCH/shepard_interp_1d.o"
