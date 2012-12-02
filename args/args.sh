#!/bin/bash
#
g++ -c -g args.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling args.cpp"
  exit
fi
rm compiler.txt
#
g++ args.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading args.o"
  exit
fi
#
rm args.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/args
#
echo "Program installed as ~/bincpp/$ARCH/args"
