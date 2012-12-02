#!/bin/bash
#
g++ -c -g triangulation_rcm.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangulation_rcm.cpp"
  exit
fi
rm compiler.txt
#
g++ triangulation_rcm.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangulation_rcm.o."
  exit
fi
#
rm triangulation_rcm.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/triangulation_rcm
#
echo "Executable installed as ~/bincpp/$ARCH/triangulation_rcm"
