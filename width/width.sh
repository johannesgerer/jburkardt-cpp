#!/bin/bash
#
g++ -c width.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling width.cpp"
  exit
fi
#
g++ width.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading width.o."
  exit
fi
#
rm width.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/width
#
echo "Executable installed as ~/bincpp/$ARCH/width"
