#!/bin/bash
#
g++ -c -g -I$HOME/include -I/opt/local/include off2plc.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling off2plc.cpp"
  exit
fi
rm compiler.txt
#
g++ off2plc.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading off2plc.o."
  exit
fi
#
rm off2plc.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/off2plc
#
echo "Executable installed as ~/bincpp/$ARCH/off2plc"
