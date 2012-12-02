#!/bin/bash
#
g++ -c -g -I/$HOME/include combo_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling combo_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ combo_prb.o /$HOME/libcpp/$ARCH/combo.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading combo_prb.o."
  exit
fi
#
rm combo_prb.o
#
mv a.out combo_prb
./combo_prb > combo_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running combo_prb."
  exit
fi
rm combo_prb
#
echo "Program output written to combo_prb_output.txt"
