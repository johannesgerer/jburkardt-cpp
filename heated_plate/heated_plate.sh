#!/bin/bash
#
g++ -c heated_plate.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling heated_plate.cpp"
  exit
fi
#
g++ heated_plate.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking heated_plate.o"
  exit
fi
#
rm heated_plate.o
mv a.out ~/bincpp/$ARCH/heated_plate
#
echo "Executable installed as ~/bincpp/$ARCH/heated_plate"
