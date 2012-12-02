#!/bin/bash
#
g++ -c -g -I/$HOME/include fsu_quality_problems.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fsu_quality_problems.cpp"
  exit
fi
rm compiler.txt
#
g++ fsu_quality_problems.o \
  /$HOME/libcpp/$ARCH/fsu_quality.o \
  /$HOME/libcpp/$ARCH/fsu_sub.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fsu_quality_problems.o."
  exit
fi
#
rm fsu_quality_problems.o
#
mv a.out fsu_quality_problems
./fsu_quality_problems > fsu_quality_problems.out
if [ $? -ne 0 ]; then
  echo "Errors running fsu_quality_problems"
  exit
fi
rm fsu_quality_problems
#
echo "The fsu_quality_problems program has been executed."
