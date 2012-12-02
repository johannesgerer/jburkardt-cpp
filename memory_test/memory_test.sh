#!/bin/bash
#
g++ -c -g -I $HOME/include memory_test.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling memory_test.cpp"
  exit
fi
rm compiler.txt
#
g++ memory_test.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading memory_test.o."
  exit
fi
#
rm memory_test.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/memory_test
#
echo "Executable installed as ~/bincpp/$ARCH/memory_test"
