#!/bin/bash
#
g++ -c -g -I/$HOME/include shepard_interp_2d_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling shepard_interp_2d_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ shepard_interp_2d_prb.o /$HOME/libcpp/$ARCH/shepard_interp_2d.o \
                            /$HOME/libcpp/$ARCH/test_interp_2d.o \
                            /$HOME/libcpp/$ARCH/r8lib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading shepard_interp_2d_prb.o."
  exit
fi
#
rm shepard_interp_2d_prb.o
#
mv a.out shepard_interp_2d_prb
./shepard_interp_2d_prb > shepard_interp_2d_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running shepard_interp_2d_prb."
  exit
fi
rm shepard_interp_2d_prb
#
echo "Program output written to shepard_interp_2d_prb_output.txt"
