#!/bin/bash
#
g++ -c -g i2_binary_to_ascii.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling i2_binary_to_ascii.cpp"
  exit
fi
rm compiler.txt
#
g++ i2_binary_to_ascii.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading i2_binary_to_ascii.o"
  exit
fi
rm i2_binary_to_ascii.o
#
mv a.out ~/bincpp/$ARCH/i2_binary_to_ascii
#
echo "Executable installed as ~/bincpp/$ARCH/i2_binary_to_ascii"
