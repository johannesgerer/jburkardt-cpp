#!/bin/bash
#
g++ -c -g -I/$HOME/include sgmg_index_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sgmg_index_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ sgmg_index_prb.o /$HOME/libcpp/$ARCH/sgmg.o /$HOME/libcpp/$ARCH/sandia_rules.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sgmg_index_prb.o."
  exit
fi
#
rm sgmg_index_prb.o
#
mv a.out sgmg_index_prb
./sgmg_index_prb > sgmg_index_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sgmg_index_prb."
  exit
fi
rm sgmg_index_prb
#
echo "Program output written to sgmg_index_prb_output.txt"
