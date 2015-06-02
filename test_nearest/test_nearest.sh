#!/bin/bash
#
g++ -c -g test_nearest.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_nearest.cpp"
  exit
fi
rm compiler.txt
#
g++ test_nearest.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading test_nearest.o"
  exit
fi
rm test_nearest.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/test_nearest
#
echo "Executable installed as ~/bincpp/$ARCH/test_nearest"
