#!/bin/bash
#
g++ -c -g problem1.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling problem1.cpp"
  exit
fi
rm compiler.txt
#
g++ problem1.o $HOME/libcpp/$ARCH/fd1d_heat_steady.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading problem1.o"
  exit
fi
rm problem1.o
#
mv a.out problem1
./problem1 > problem1_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running problem1"
  exit
fi
rm problem1
#
echo "Test program output written to problem1_output.txt."
