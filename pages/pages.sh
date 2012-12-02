#!/bin/bash
#
g++ -c -g pages.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pages.cpp"
  exit
fi
rm compiler.txt
#
g++ pages.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pages.o."
  exit
fi
#
rm pages.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/pages
#
echo "Executable installed as ~/bincpp/$ARCH/pages"
