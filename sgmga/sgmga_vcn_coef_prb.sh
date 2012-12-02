#!/bin/bash
#
g++ -c -g -I/$HOME/include sgmga_vcn_coef_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sgmga_vcn_coef_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ sgmga_vcn_coef_prb.o /$HOME/libcpp/$ARCH/sgmga.o /$HOME/libcpp/$ARCH/sandia_rules.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sgmga_vcn_coef_prb.o."
  exit
fi
#
rm sgmga_vcn_coef_prb.o
#
mv a.out sgmga_vcn_coef_prb
./sgmga_vcn_coef_prb > sgmga_vcn_coef_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sgmga_vcn_coef_prb."
  exit
fi
rm sgmga_vcn_coef_prb
#
echo "Program output written to sgmga_vcn_coef_prb_output.txt"
