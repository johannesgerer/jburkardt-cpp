#!/bin/bash
#
g++ -c decomment.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling decomment.cpp"
  exit
fi
rm compiler.txt
#
g++ decomment.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading decomment.o."
  exit
fi
#
rm decomment.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/decomment
#
echo "Executable installed as ~/bincpp/$ARCH/decomment"
