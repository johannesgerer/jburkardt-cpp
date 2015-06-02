#!/bin/bash
#
g++ -c -g -I/$HOME/include problem2.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling problem2.cpp"
  exit
fi
rm compiler.txt
#
g++ /$HOME/libcpp/$ARCH/dream.o \
  problem2.o \
  /$HOME/libcpp/$ARCH/pdflib.o \
  /$HOME/libcpp/$ARCH/rnglib.o -lm
#
if [ $? -ne 0 ]; then
  echo "Errors linking and loading dream.o + problem2.o + pdflib.o + rnglib.o"
  exit
fi
#
rm problem2.o
#
mv a.out problem2
./problem2 > problem2_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running problem2."
  exit
fi
rm problem2
#
echo "Program output written to problem2_output.txt"
