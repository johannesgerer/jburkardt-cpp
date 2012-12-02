#!/bin/bash
#
g++ -c -g -I/$HOME/include piecewise_linear_product_integral_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling piecewise_linear_product_integral_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ piecewise_linear_product_integral_prb.o /$HOME/libcpp/$ARCH/piecewise_linear_product_integral.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading piecewise_linear_product_integral_prb.o."
  exit
fi
#
rm piecewise_linear_product_integral_prb.o
#
mv a.out piecewise_linear_product_integral_prb
./piecewise_linear_product_integral_prb > piecewise_linear_product_integral_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running piecewise_linear_product_integral_prb."
  exit
fi
rm piecewise_linear_product_integral_prb
#
echo "Program output written to piecewise_linear_product_integral_prb_output.txt"
