#!/bin/bash
#
g++ -c -I$HOME/include ppma_to_ppmb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ppma_to_ppmb.cpp"
  exit
fi
rm compiler.txt
#
g++ ppma_to_ppmb.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ppma_to_ppmb.o"
  exit
fi
#
rm ppma_to_ppmb.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/ppma_to_ppmb
#
echo "A new version of ppma_to_ppmb has been created."
