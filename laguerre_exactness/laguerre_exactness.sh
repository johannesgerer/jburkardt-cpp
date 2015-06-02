#!/bin/bash
#
g++ -c laguerre_exactness.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling laguerre_exactness.cpp"
  exit
fi
rm compiler.txt
#
g++ laguerre_exactness.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading laguerre_exactness.o"
  exit
fi
rm laguerre_exactness.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/laguerre_exactness
#
echo "Executable installed as ~/bincpp/$ARCH/laguerre_exactness"
