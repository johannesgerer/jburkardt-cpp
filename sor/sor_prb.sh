#!/bin/bash
#
g++ -c -g -I/$HOME/include sor_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sor_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ sor_prb.o /$HOME/libcpp/$ARCH/sor.o ~/libc/$ARCH/gnuplot_i.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sor_prb.o."
  exit
fi
#
rm sor_prb.o
#
mv a.out sor_prb
./sor_prb > sor_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sor_prb."
  exit
fi
rm sor_prb
#
echo "Program output written to sor_prb_output.txt"
