#!/bin/bash
#
g++ -c -g  -I/$HOME/include simple_rkf45.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling simple_rkf45.cpp"
  exit
fi
rm compiler.txt
#
g++ simple_rkf45.o /$HOME/libcpp/$ARCH/rkf45.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading simple_rkf45.o"
  exit
fi
rm simple_rkf45.o
#
mv a.out simple_rkf45
./simple_rkf45 > simple_rkf45_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running simple_rkf45"
  exit
fi
rm simple_rkf45
#
echo "Test program output written to simple_rkf45_output.txt."
