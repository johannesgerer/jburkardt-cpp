#!/bin/bash
#
g++ -c -g -I$HOME/include fsu_cvt_standalone.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fsu_cvt_standalone.cpp"
  exit
fi
rm compiler.txt
#
g++ fsu_cvt_standalone.o \
  $HOME/libcpp/$ARCH/fsu_cvt.o \
  $HOME/libcpp/$ARCH/fsu_halton.o \
  $HOME/libcpp/$ARCH/fsu_sub.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fsu_cvt_standalone.o"
  echo "                         + fsu_cvt.o."
  echo "                         + fsu_halton.o."
  echo "                         + fsu_sub.o."
  exit
fi
#
rm fsu_cvt_standalone.o
mv a.out ~/bincpp/$ARCH/fsu_cvt_standalone
#
echo "A new version of fsu_cvt_standalone has been created."
