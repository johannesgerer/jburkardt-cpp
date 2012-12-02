#!/bin/bash
#
g++ -c -I$HOME/include pgma_to_pgmb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pgma_to_pgmb.cpp"
  exit
fi
rm compiler.txt
#
g++ pgma_to_pgmb.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pgma_to_pgmb.o"
  exit
fi
#
rm pgma_to_pgmb.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/pgma_to_pgmb
#
echo "Executable installed as ~/bincpp/$ARCH/pgma_to_pgmb"
