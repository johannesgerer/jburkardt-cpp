#!/bin/bash
#
g++ -c fd1d_advection_lax_wendroff.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling fd1d_advection_lax_wendroff.cpp"
  exit
fi
#
g++ fd1d_advection_lax_wendroff.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking fd1d_advection_lax_wendroff.o"
  exit
fi
#
rm fd1d_advection_lax_wendroff.o
mv a.out ~/bincpp/$ARCH/fd1d_advection_lax_wendroff
#
echo "Executable installed as ~/bincpp/$ARCH/fd1d_advection_lax_wendroff"
