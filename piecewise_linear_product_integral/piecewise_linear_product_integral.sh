#!/bin/bash
#
cp piecewise_linear_product_integral.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include piecewise_linear_product_integral.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling piecewise_linear_product_integral.cpp"
  exit
fi
rm compiler.txt
#
mv piecewise_linear_product_integral.o ~/libcpp/$ARCH/piecewise_linear_product_integral.o
#
echo "Library installed as ~/libcpp/$ARCH/piecewise_linear_product_integral.o"
