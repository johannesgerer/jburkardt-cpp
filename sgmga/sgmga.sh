#!/bin/bash
#
cp sgmga.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include sgmga.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sgmga.cpp"
  exit
fi
rm compiler.txt
#
mv sgmga.o ~/libcpp/$ARCH/sgmga.o
#
echo "Library installed as ~/libcpp/$ARCH/sgmga.o"
