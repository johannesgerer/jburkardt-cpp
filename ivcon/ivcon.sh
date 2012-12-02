#!/bin/bash
#
g++ -c ivcon.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ivcon.cpp"
  exit
fi
rm compiler.txt
#
g++ ivcon.o -lm
if [ $? -ne 0 ]; then
  ercho "Errors linking and loading ivcon.o"
  exit
fi
#
rm ivcon.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/ivcon
#
echo "Executable installed as ~/bincpp/$ARCH/ivcon"
