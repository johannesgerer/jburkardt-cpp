#!/bin/bash
#
g++ -c catalog.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling catalog.cpp"
  exit
fi
rm compiler.txt
#
g++ catalog.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading catalog.o."
  exit
fi
#
rm catalog.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/catalog
#
echo "Executable installed as ~/bincpp/$ARCH/catalog"
