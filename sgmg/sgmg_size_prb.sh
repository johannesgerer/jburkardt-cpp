#!/bin/bash
#
g++ -c -g -I/$HOME/include sgmg_size_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sgmg_size_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ sgmg_size_prb.o /$HOME/libcpp/$ARCH/sgmg.o /$HOME/libcpp/$ARCH/sandia_rules.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sgmg_size_prb.o."
  exit
fi
#
rm sgmg_size_prb.o
#
mv a.out sgmg_size_prb
./sgmg_size_prb > sgmg_size_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sgmg_size_prb."
  exit
fi
rm sgmg_size_prb
#
echo "Program output written to sgmg_size_prb_output.txt"
