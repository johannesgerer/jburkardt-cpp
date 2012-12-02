#!/bin/bash
#
g++ -c -g -I/$HOME/include sgmg_point_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sgmg_point_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ sgmg_point_prb.o /$HOME/libcpp/$ARCH/sgmg.o /$HOME/libcpp/$ARCH/sandia_rules.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sgmg_point_prb.o."
  exit
fi
#
rm sgmg_point_prb.o
#
mv a.out sgmg_point_prb
./sgmg_point_prb > sgmg_point_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sgmg_point_prb."
  exit
fi
rm sgmg_point_prb
#
echo "Program output written to sgmg_point_prb_output.txt"
