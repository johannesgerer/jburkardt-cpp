#!/bin/bash
#
g++ -c -g -I/$HOME/include latin_center_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling latin_center_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ latin_center_prb.o /$HOME/libcpp/$ARCH/latin_center.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading latin_center_prb.o."
  exit
fi
#
rm latin_center_prb.o
#
mv a.out latin_center_prb
./latin_center_prb > latin_center_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running latin_center_prb."
  exit
fi
rm latin_center_prb
#
echo "Program output written to latin_center_prb_output.txt"
