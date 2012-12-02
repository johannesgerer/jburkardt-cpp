#!/bin/bash
#
cp partition_problem.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include partition_problem.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling partition_problem.cpp"
  exit
fi
rm compiler.txt
#
mv partition_problem.o ~/libcpp/$ARCH/partition_problem.o
#
echo "Library installed as ~/libcpp/$ARCH/partition_problem.o"
