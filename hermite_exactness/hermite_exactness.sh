#!/bin/bash
#
g++ -c hermite_exactness.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling hermite_exactness.cpp"
  exit
fi
#
g++ hermite_exactness.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading hermite_exactness.o"
  exit
fi
rm hermite_exactness.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/hermite_exactness
#
echo "Executable installed as ~/bincpp/$ARCH/hermite_exactness"
