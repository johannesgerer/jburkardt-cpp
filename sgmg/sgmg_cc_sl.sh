#!/bin/bash
#
g++ -c -g -I/$HOME/include sgmg_cc_sl.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sgmg_cc_sl.cpp"
  exit
fi
rm compiler.txt
#
g++ sgmg_cc_sl.o /$HOME/libcpp/$ARCH/sgmg.o /$HOME/libcpp/$ARCH/sandia_rules.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sgmg_cc_sl.o."
  exit
fi
#
rm sgmg_cc_sl.o
#
mv a.out sgmg_cc_sl
./sgmg_cc_sl > sgmg_cc_sl_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sgmg_cc_sl."
  exit
fi
rm sgmg_cc_sl
#
echo "Program output written to sgmg_cc_sl_output.txt"

