#!/bin/bash
#
g++ -c -g -I$HOME/include fsu_quality_standalone.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fsu_quality_standalone.cpp"
  exit
fi
rm compiler.txt
#
g++ fsu_quality_standalone.o \
  $HOME/libcpp/$ARCH/fsu_quality.o \
  $HOME/libcpp/$ARCH/fsu_sub.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fsu_quality_standalone.o"
  echo "                         + fsu_quality.o."
  echo "                         + fsu_sub.o."
  exit
fi
#
rm fsu_quality_standalone.o
mv a.out ~/bincpp/$ARCH/fsu_quality_standalone
#
echo "A new version of fsu_quality_standalone has been created."
