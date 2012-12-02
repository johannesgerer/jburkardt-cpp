#!/bin/bash
#
g++ -c -g -I/$HOME/include fsu_latinize_problems.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fsu_latinize_problems.cpp"
  exit
fi
rm compiler.txt
#
g++ fsu_latinize_problems.o \
  /$HOME/libcpp/$ARCH/fsu_latinize.o \
  /$HOME/libcpp/$ARCH/fsu_halton.o \
  /$HOME/libcpp/$ARCH/fsu_sub.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fsu_latinize_problems.o."
  exit
fi
#
rm fsu_latinize_problems.o
#
mv a.out fsu_latinize_problems
./fsu_latinize_problems > fsu_latinize_problems.out
if [ $? -ne 0 ]; then
  echo "Errors running fsu_latinize_problems"
  exit
fi
rm fsu_latinize_problems
#
echo "The fsu_latinize_problems program has been executed."
