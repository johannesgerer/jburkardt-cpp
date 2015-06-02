#!/bin/bash
#
g++ -c -g life_serial.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling life_serial.cpp"
  exit
fi
rm compiler.txt
#
g++ life_serial.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading life_serial.o."
  exit
fi
#
rm life_serial.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/life_serial
#
echo "Executable installed as ~/bincpp/$ARCH/life_serial"
