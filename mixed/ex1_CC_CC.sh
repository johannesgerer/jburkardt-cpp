#!/bin/bash
#
#  Compile the C++ file.
#
CC -c ex1_main.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ex1_main.cpp"
  exit
fi
rm compiler.txt
#
#  Compile the C file.
#
CC -c ex1_sub.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ex1_sub.c"
  exit
fi
rm compiler.txt
#
#  Use the CC program to load.
#  In some cases, the C libraries may need to be included,
#  using an argument like "-lc", as well as the C math libraries,
#  using an argument like "-lm".
#
CC ex1_main.o ex1_sub.o -lc -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ex1_main.o + ex1_sub.o"
  exit
fi
#
mv a.out ex1_CC_CC
rm *.o
#
#  Run the program.
#
ex1_CC_CC > ex1_CC_CC_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ex1_CC_CC."
  exit
fi
rm ex1_CC_CC
#
echo "Program output written to ex1_CC_CC_output.txt"
