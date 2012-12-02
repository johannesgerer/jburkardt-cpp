#!/bin/bash
#
g++ -c -g gfplys.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling gfplys.cpp"
  exit
fi
rm compiler.txt
#
g++ gfplys.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading gfplys.o."
  exit
fi
#
rm gfplys.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/gfplys
#
echo "Executable installed as ~/bincpp/$ARCH/gfplys"
