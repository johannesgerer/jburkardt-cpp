#!/bin/bash
#
g++ -c -g fem2d_poisson_rectangle_linear.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem2d_poisson_rectangle_linear.cpp"
  exit
fi
rm compiler.txt
#
g++ fem2d_poisson_rectangle_linear.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem2d_poisson_rectangle_linear.o"
  exit
fi
#
rm fem2d_poisson_rectangle_linear.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/fem2d_poisson_rectangle_linear
#
echo "Executable installed as ~/bincpp/$ARCH/fem2d_poisson_rectangle_linear"
