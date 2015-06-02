#!/bin/bash
#
g++ -c fd1d_advection_lax.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fd1d_advection_lax.cpp"
  exit
fi
rm compiler.txt
#
g++ fd1d_advection_lax.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking fd1d_advection_lax.o"
  exit
fi
#
rm fd1d_advection_lax.o
mv a.out ~/bincpp/$ARCH/fd1d_advection_lax
#
echo "Executable installed as ~/bincpp/$ARCH/fd1d_advection_lax"
