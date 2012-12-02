#!/bin/bash
#
g++ -c -g -I /$HOME/include calendar_rd.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling calendar_rd.cpp"
  exit
fi
rm compiler.txt
#
g++ calendar_rd.o -lm 
if [ $? -ne 0 ]; then
  echo "Errors linking and loading calendar_rd.o."
  exit
fi
rm calendar_rd.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/calendar_rd
#
echo "Executable installed as ~/bincpp/$ARCH/calendar_rd"
