#!/bin/bash
#
g++ -c fd1d_advection_diffusion_steady.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fd1d_advection_diffusion_steady.cpp"
  exit
fi
rm compiler.txt
#
g++ fd1d_advection_diffusion_steady.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking fd1d_advection_diffusion_steady.o"
  exit
fi
#
rm fd1d_advection_diffusion_steady.o
mv a.out ~/bincpp/$ARCH/fd1d_advection_diffusion_steady
#
echo "Executable installed as ~/bincpp/$ARCH/fd1d_advection_diffusion_steady"
