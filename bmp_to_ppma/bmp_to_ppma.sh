#!/bin/bash
#
g++ -c -I$HOME/include bmp_to_ppma.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling bmp_to_ppma.cpp"
  exit
fi
rm compiler.txt
#
g++ bmp_to_ppma.o ~/libcpp/$ARCH/bmp_io.o ~/libcpp/$ARCH/ppma_io.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading bmp_to_ppma.o"
  exit
fi
#
rm bmp_to_ppma.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/bmp_to_ppma
#
echo "Executable installed as ~/bincpp/$ARCH/bmp_to_ppma"
