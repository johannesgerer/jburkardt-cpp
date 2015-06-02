#!/bin/bash
#
g++ -c -I$HOME/include hypercube_exactness.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling hypercube_exactness.cpp"
  exit
fi
#
g++ hypercube_exactness.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading hypercube_exactness.o."
  exit
fi
#
rm hypercube_exactness.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/hypercube_exactness
#
echo "Executable installed as ~/bincpp/$ARCH/hypercube_exactness"
