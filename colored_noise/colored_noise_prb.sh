#!/bin/bash
#
g++ -c -g -I/$HOME/include colored_noise_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling colored_noise_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ colored_noise_prb.o /$HOME/libcpp/$ARCH/colored_noise.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading colored_noise_prb.o."
  exit
fi
#
rm colored_noise_prb.o
#
mv a.out colored_noise_prb
./colored_noise_prb > colored_noise_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running colored_noise_prb."
  exit
fi
rm colored_noise_prb
#
echo "Program output written to colored_noise_prb_output.txt"
