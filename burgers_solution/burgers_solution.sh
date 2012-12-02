#!/bin/bash
#
cp burgers_solution.hpp /$HOME/include
#
g++ -c -g burgers_solution.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling burgers_solution.cpp"
  exit
fi
rm compiler.txt
#
mv burgers_solution.o ~/libcpp/$ARCH/burgers_solution.o
#
echo "Library installed as ~/libcpp/$ARCH/burgers_solution.o"
