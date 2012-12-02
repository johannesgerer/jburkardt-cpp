#!/bin/bash
#
cp task_division.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include task_division.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling task_division.cpp"
  exit
fi
rm compiler.txt
#
mv task_division.o ~/libcpp/$ARCH/task_division.o
#
echo "Library installed as ~/libcpp/$ARCH/task_division.o"
