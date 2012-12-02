#!/bin/bash
#
g++ -c -g -I/$HOME/include sgmga_size_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sgmga_size_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ sgmga_size_prb.o /$HOME/libcpp/$ARCH/sgmga.o /$HOME/libcpp/$ARCH/sandia_rules.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sgmga_size_prb.o."
  exit
fi
#
rm sgmga_size_prb.o
#
mv a.out sgmga_size_prb
./sgmga_size_prb > sgmga_size_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sgmga_size_prb."
  exit
fi
rm sgmga_size_prb
#
echo "Program output written to sgmga_size_prb_output.txt"
