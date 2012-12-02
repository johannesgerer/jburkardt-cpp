#!/bin/bash
#
g++ -c reword.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling reword.cpp"
  exit
fi
rm compiler.txt
#
g++ reword.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading reword.o."
  exit
fi
#
rm reword.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/reword
#
echo "Executable installed as ~/bincpp/$ARCH/reword"
