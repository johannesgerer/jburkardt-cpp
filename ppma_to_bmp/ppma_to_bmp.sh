#!/bin/bash
#
g++ -c -I$HOME/include ppma_to_bmp.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ppma_to_bmp.cpp"
  exit
fi
rm compiler.txt
#
g++ ppma_to_bmp.o ~/libcpp/$ARCH/bmp_io.o ~/libcpp/$ARCH/ppma_io.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ppma_to_bmp.o"
  exit
fi
#
rm ppma_to_bmp.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/ppma_to_bmp
#
echo "Executable installed as ~/bincpp/$ARCH/ppma_to_bmp"
