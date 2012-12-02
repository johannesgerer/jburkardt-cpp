#!/bin/bash
#
g++ -c recomment.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling recomment.cpp"
  exit
fi
rm compiler.txt
#
g++ recomment.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading recomment.o."
  exit
fi
#
rm recomment.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/recomment
#
echo "Executable installed as ~/bin/$ARCH/recomment"
