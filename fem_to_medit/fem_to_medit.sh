#!/bin/bash
#
g++ -c fem_to_medit.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling fem_to_medit.cpp"
  exit
fi
#
g++ fem_to_medit.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem_to_medit.o"
  exit
fi
#
rm fem_to_medit.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/fem_to_medit
#
echo "Program installed as ~/bincpp/$ARCH/fem_to_medit"
