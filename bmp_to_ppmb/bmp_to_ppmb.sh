#!/bin/bash
#
g++ -c -I $HOME/include bmp_to_ppmb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling bmp_to_ppmb.cpp"
  exit
fi
rm compiler.txt
#
g++ bmp_to_ppmb.o ~/libcpp/$ARCH/bmp_io.o ~/libcpp/$ARCH/ppmb_io.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading bmp_to_ppmb.o"
  exit
fi
#
rm bmp_to_ppmb.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/bmp_to_ppmb
#
echo "Executable installed as ~/bincpp/$ARCH/bmp_to_ppmb"
