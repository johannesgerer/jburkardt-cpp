#!/bin/bash
#
g++ -c rot13.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling rot13.cpp"
  exit
fi
rm compiler.txt
#
g++ rot13.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading rot13.o."
  exit
fi
#
rm rot13.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/rot13
#
echo "The executable was installed as ~/bincpp/$ARCH/rot13"
