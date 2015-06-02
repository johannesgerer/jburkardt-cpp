#!/bin/bash
#
g++ -c -g allocatable.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling allocatable.cpp"
  exit
fi
rm compiler.txt
#
g++ allocatable.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading allocatable.o"
  exit
fi
rm allocatable.o
#
mv a.out allocatable
./allocatable > allocatable_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running allocatable"
  exit
fi
rm allocatable
#
echo "Program output written to allocatable_output.txt"
