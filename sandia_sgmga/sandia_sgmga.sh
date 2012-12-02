#!/bin/bash
#
cp sandia_sgmga.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include sandia_sgmga.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sandia_sgmga.cpp"
  exit
fi
rm compiler.txt
#
mv sandia_sgmga.o ~/libcpp/$ARCH/sandia_sgmga.o
#
echo "Library installed as ~/libcpp/$ARCH/sandia_sgmga.o"
