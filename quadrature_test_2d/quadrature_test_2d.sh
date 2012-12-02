#!/bin/bash
#
g++ -c -g -I/$HOME/include quadrature_test_2d.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling quadrature_test_2d.cpp"
  exit
fi
rm compiler.txt
#
g++ quadrature_test_2d.o /$HOME/libcpp/$ARCH/test_int_2d.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading quadrature_test_2d.o + testnint.o"
  exit
fi
#
rm quadrature_test_2d.o
#
mv a.out $HOME/bincpp/$ARCH/quadrature_test_2d
#
echo "Executable installed as ~/bincpp/$ARCH/quadrature_test_2d"
