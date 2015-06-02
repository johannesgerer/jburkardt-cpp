#!/bin/bash
#
g++ -c -I/$HOME/include halton_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling halton_prb.cpp"
  exit
fi
#
g++ halton_prb.o /$HOME/libcpp/$ARCH/halton.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading halton_prb.o."
  exit
fi
#
rm halton_prb.o
#
mv a.out halton_prb
./halton_prb > halton_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running halton_prb."
  exit
fi
rm halton_prb
#
echo "Program output written to halton_prb_output.txt"
