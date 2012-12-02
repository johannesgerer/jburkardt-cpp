#!/bin/bash
#
g++ -c uncontrol.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling uncontrol.cpp"
  exit
fi
rm compiler.txt
#
g++ uncontrol.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading uncontrol.o."
  exit
fi
#
rm uncontrol.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/uncontrol
#
echo "Executable installed as ~/bincpp/$ARCH/uncontrol"
