#!/bin/bash
#
#  Compile the C++ file.
#
g++ -c ex1_main.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ex1_main.cpp"
  exit
fi
rm compiler.txt
#
#  Compile the C file.
#
gcc -c ex1_sub.c
if [ $? -ne 0 ]; then
  echo "Errors compiling ex1_sub.c"
  exit
fi
#
#  Use the g++ program to load.
#  In some cases, the C libraries may need to be included,
#  using an argument like "-lc", as well as the C math libraries,
#  using an argument like "-lm".
#
g++ ex1_main.o ex1_sub.o -lc -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ex1_main.o + ex1_sub.o"
  exit
fi
#
mv a.out ex1_g++_gcc
rm *.o
#
#  Run the program.
#
ex1_g++_gcc > ex1_g++_gcc_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ex1_g++_gcc."
  exit
fi
rm ex1_g++_gcc
#
echo "Program output written to ex1_g++_gcc_output.txt"
