#!/bin/bash
#
g++ -c limits.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling limits.cpp"
  exit
fi
rm compiler.txt
#
g++ limits.o
if [ $? -ne 0 ]; then
  echo "Errors loading limits.o"
  exit
fi
rm limits.o
#
mv a.out limits
./limits > limits_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running limits"
  exit
fi
rm limits
#
echo "Program output written to limits_output.txt"
