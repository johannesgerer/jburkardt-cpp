#!/bin/bash
#
g++ -c -g -I/$HOME/include problem1_main.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling problem1_main.cpp"
  exit
fi
rm compiler.txt
#
g++ -c -g -I/$HOME/include problem1.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling problem1.cpp"
  exit
fi
rm compiler.txt
#
g++ problem1_main.o \
    problem1.o \
  /$HOME/libcpp/$ARCH/pdflib.o \
  /$HOME/libcpp/$ARCH/rnglib.o \
  -lm
if [ $? -ne 0 ]; then
  echo "Errors loading problem1_main.o + problem1.o + pdflib.o + rnglib.o"
  exit
fi
#
rm problem1_main.o
rm problem1.o
#
mv a.out problem1_main
./problem1_main > problem1_main_output.txt
rm problem1_main
#
echo "Program output written to problem1_main_output.txt"
