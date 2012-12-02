#!/bin/bash
#
g++ -c png_crc.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling png_crc.cpp"
  exit
fi
rm compiler.txt
#
g++ png_crc.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading png_crc.o."
  exit
fi
#
rm png_crc.o
mv a.out ~/bincpp/$ARCH/png_crc
#
echo "Executable installed as ~/bincpp/$ARCH/png_crc"
