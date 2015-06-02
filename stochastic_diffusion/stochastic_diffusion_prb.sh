#!/bin/bash
#
g++ -c -g -I/$HOME/include stochastic_diffusion_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling stochastic_diffusion_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ stochastic_diffusion_prb.o /$HOME/libcpp/$ARCH/stochastic_diffusion.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading stochastic_diffusion_prb.o."
  exit
fi
#
rm stochastic_diffusion_prb.o
#
mv a.out stochastic_diffusion_prb
./stochastic_diffusion_prb > stochastic_diffusion_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running stochastic_diffusion_prb."
  exit
fi
rm stochastic_diffusion_prb
#
echo "Program output written to stochastic_diffusion_prb_output.txt"
