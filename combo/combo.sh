#!/bin/bash
#
cp combo.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include combo.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling combo.cpp"
  exit
fi
rm compiler.txt
#
mv combo.o ~/libcpp/$ARCH/combo.o
#
echo "Library installed as ~/libcpp/$ARCH/combo.o"
