#!/bin/bash
#
g++ -c -g ccvt_box.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ccvt_box.cpp"
  exit
fi
rm compiler.txt
#
g++ ccvt_box.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ccvt_box.o"
  exit
fi
rm ccvt_box.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/ccvt_box
#
echo "Executable installed as ~/bincpp/$ARCH/ccvt_box"
