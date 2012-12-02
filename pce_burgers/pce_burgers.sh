#!/bin/bash
#
g++ -c pce_burgers.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pce_burgers.cpp"
  exit
fi
rm compiler.txt
#
g++ pce_burgers.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pce_burgers.o"
  exit
fi
rm pce_burgers.o
#
mv a.out ~/bincpp/$ARCH/pce_burgers
#
echo "Executable installed as ~/bincpp/$ARCH/pce_burgers"
