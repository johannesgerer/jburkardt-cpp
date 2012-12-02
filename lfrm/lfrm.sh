#!/bin/bash
#
g++ -c lfrm.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "There were errors while compiling lfrm.cpp"
  exit
fi
rm compiler.txt
#
g++ lfrm.o
if [ $? -ne 0 ]; then
  exit
fi
#
rm lfrm.o
mv a.out ~/bincpp/$ARCH/lfrm
#
echo "Executable installed as ~/bincpp/$ARCH/lfrm."
