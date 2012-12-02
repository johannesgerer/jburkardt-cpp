#!/bin/bash
#
g++ -c -g inout.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling inout.cpp"
  exit
fi
rm compiler.txt
#
g++ ~/libcpp/$ARCH/free_fem_stokes.o inout.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading free_fem_stokes.o + inout.o"
  exit
fi
rm inout.o
#
chmod ugo+x a.out
mv a.out inout
./inout nodes6.txt triangles6.txt > inout_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running inout."
  exit
fi
rm inout
#
echo "The inout program has been executed."
