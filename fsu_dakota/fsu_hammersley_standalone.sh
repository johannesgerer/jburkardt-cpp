#!/bin/bash
#
g++ -c -g -I$HOME/include fsu_hammersley_standalone.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fsu_hammersley_standalone.cpp"
  exit
fi
rm compiler.txt
#
g++ fsu_hammersley_standalone.o \
  $HOME/libcpp/$ARCH/fsu_hammersley.o \
  $HOME/libcpp/$ARCH/fsu_sub.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fsu_hammersley_standalone.o"
  echo "                         + fsu_hammersley.o."
  echo "                         + fsu_sub.o."
  exit
fi
#
rm fsu_hammersley_standalone.o
mv a.out ~/bincpp/$ARCH/fsu_hammersley_standalone
#
echo "A new version of fsu_hammersley_standalone has been created."
