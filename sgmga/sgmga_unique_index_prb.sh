#!/bin/bash
#
g++ -c -g -I/$HOME/include sgmga_unique_index_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sgmga_unique_index_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ sgmga_unique_index_prb.o /$HOME/libcpp/$ARCH/sgmga.o /$HOME/libcpp/$ARCH/sandia_rules.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sgmga_unique_index_prb.o."
  exit
fi
#
rm sgmga_unique_index_prb.o
#
mv a.out sgmga_unique_index_prb
./sgmga_unique_index_prb > sgmga_unique_index_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sgmga_unique_index_prb."
  exit
fi
rm sgmga_unique_index_prb
#
echo "Program output written to sgmga_unique_index_prb_output.txt"
