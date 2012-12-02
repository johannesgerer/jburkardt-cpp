#!/bin/bash
#
g++ -c cr2crlf.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cr2crlf.cpp"
  exit
fi
rm compiler.txt
#
g++ cr2crlf.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cr2crlf.o."
  exit
fi
#
rm cr2crlf.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/cr2crlf
#
echo "Executable installed as ~/bincpp/$ARCH/cr2crlf."
