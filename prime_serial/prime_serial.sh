#!/bin/bash
#
cp prime_serial.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include prime_serial.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling prime_serial.cpp"
  exit
fi
rm compiler.txt
#
mv prime_serial.o ~/libcpp/$ARCH/prime_serial.o
#
echo "Library installed as ~/libcpp/$ARCH/prime_serial.o"
