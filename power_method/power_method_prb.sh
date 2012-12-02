#!/bin/bash
#
g++ -c -g -I/$HOME/include power_method_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling power_method_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ power_method_prb.o /$HOME/libcpp/$ARCH/power_method.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading power_method_prb.o."
  exit
fi
#
rm power_method_prb.o
#
mv a.out power_method_prb
./power_method_prb > power_method_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running power_method_prb."
  exit
fi
rm power_method_prb
#
echo "Program output written to power_method_prb_output.txt"
