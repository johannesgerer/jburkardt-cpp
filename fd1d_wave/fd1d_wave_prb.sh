#!/bin/bash
#
g++ -c -g -I/$HOME/include fd1d_wave_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fd1d_wave_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ fd1d_wave_prb.o /$HOME/libcpp/$ARCH/fd1d_wave.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fd1d_wave_prb.o."
  exit
fi
#
rm fd1d_wave_prb.o
#
mv a.out fd1d_wave_prb
./fd1d_wave_prb > fd1d_wave_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running fd1d_wave_prb."
  exit
fi
rm fd1d_wave_prb
#
echo "Program output written to fd1d_wave_prb_output.txt"
