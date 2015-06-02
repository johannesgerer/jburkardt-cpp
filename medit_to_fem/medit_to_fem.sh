#!/bin/bash
#
g++ -c medit_to_fem.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling medit_to_fem.cpp"
  exit
fi
#
g++ medit_to_fem.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading medit_to_fem.o."
  exit
fi
#
rm medit_to_fem.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/medit_to_fem
#
echo "Executable installed as ~/bincpp/$ARCH/medit_to_fem"
