#!/bin/bash
#
g++ -c -g int_exactness_jacobi.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling int_exactness_jacobi.cpp"
  exit
fi
rm compiler.txt
#
g++ int_exactness_jacobi.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading int_exactness_jacobi.o"
  exit
fi
rm int_exactness_jacobi.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/int_exactness_jacobi
#
echo "Executable installed as ~/bincpp/$ARCH/int_exactness_jacobi"
