#!/bin/bash
#
g++ -c -g -I/$HOME/include sgmga_aniso_normalize_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sgmga_aniso_normalize_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ sgmga_aniso_normalize_prb.o /$HOME/libcpp/$ARCH/sgmga.o /$HOME/libcpp/$ARCH/sandia_rules.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sgmga_aniso_normalize_prb.o."
  exit
fi
#
rm sgmga_aniso_normalize_prb.o
#
mv a.out sgmga_aniso_normalize_prb
./sgmga_aniso_normalize_prb > sgmga_aniso_normalize_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sgmga_aniso_normalize_prb."
  exit
fi
rm sgmga_aniso_normalize_prb
#
echo "Program output written to sgmga_aniso_normalize_prb_output.txt"
