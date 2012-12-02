#!/bin/bash
#
g++ -c -g -I/$HOME/include fsu_cvt_problems.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fsu_cvt_problems.cpp"
  exit
fi
rm compiler.txt
#
g++ fsu_cvt_problems.o \
  /$HOME/libcpp/$ARCH/fsu_cvt.o \
  /$HOME/libcpp/$ARCH/fsu_halton.o \
  /$HOME/libcpp/$ARCH/fsu_sub.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fsu_cvt_problems.o."
  exit
fi
#
rm fsu_cvt_problems.o
#
mv a.out fsu_cvt_problems
./fsu_cvt_problems > fsu_cvt_problems.out
if [ $? -ne 0 ]; then
  echo "Errors running fsu_cvt_problems"
  exit
fi
rm fsu_cvt_problems
#
echo "The fsu_cvt_problems program has been executed."
