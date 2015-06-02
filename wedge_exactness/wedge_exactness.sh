#!/bin/bash
#
g++ -c wedge_exactness.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling wedge_exactness.cpp"
  exit
fi
#
g++ wedge_exactness.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading wedge_exactness.o"
  exit
fi
rm wedge_exactness.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/wedge_exactness
#
echo "Executable installed as ~/bincpp/$ARCH/wedge_exactness"
