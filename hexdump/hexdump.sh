#!/bin/bash
#
g++ -c hexdump.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling hexdump.cpp"
  exit
fi
#
g++ hexdump.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading hexdump.o."
  exit
fi
#
rm hexdump.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/hexdump
#
echo "Executable installed as ~/bincpp/$ARCH/hexdump"
