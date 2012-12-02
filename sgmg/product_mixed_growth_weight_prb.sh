#!/bin/bash
#
g++ -c -g -I/$HOME/include product_mixed_growth_weight_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling product_mixed_growth_weight_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ product_mixed_growth_weight_prb.o /$HOME/libcpp/$ARCH/sgmg.o /$HOME/libcpp/$ARCH/sandia_rules.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading product_mixed_growth_weight_prb.o."
  exit
fi
#
rm product_mixed_growth_weight_prb.o
#
mv a.out product_mixed_growth_weight_prb
./product_mixed_growth_weight_prb > product_mixed_growth_weight_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running product_mixed_growth_weight_prb."
  exit
fi
rm product_mixed_growth_weight_prb
#
echo "Program output written to product_mixed_growth_weight_prb_output.txt"
