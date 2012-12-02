#!/bin/bash
#
cp sgmg.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include sgmg.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sgmg.cpp"
  exit
fi
rm compiler.txt
#
mv sgmg.o ~/libcpp/$ARCH/sgmg.o
#
echo "Library installed as ~/libcpp/$ARCH/sgmg.o"
