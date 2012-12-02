#!/bin/bash
#
g++ -c -g -I$HOME/include fsu_latinize_standalone.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fsu_latinize_standalone.cpp"
  exit
fi
rm compiler.txt
#
g++ fsu_latinize_standalone.o \
  $HOME/libcpp/$ARCH/fsu_latinize.o \
  $HOME/libcpp/$ARCH/fsu_sub.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fsu_latinize_standalone.o"
  echo "                         + fsu_latinize.o."
  echo "                         + fsu_sub.o."
  exit
fi
#
rm fsu_latinize_standalone.o
mv a.out ~/bincpp/$ARCH/fsu_latinize_standalone
#
echo "A new version of fsu_latinize_standalone has been created."
