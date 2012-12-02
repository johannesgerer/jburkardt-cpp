#!/bin/bash
#
cp sandia_sgmgg.hpp $HOME/include
#
g++ -c -g -I $HOME/include sandia_sgmgg.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sandia_sgmgg.cpp"
  exit
fi
rm compiler.txt
#
mv sandia_sgmgg.o ~/libcpp/$ARCH/sandia_sgmgg.o
#
echo "Library installed as ~/libcpp/$ARCH/sandia_sgmgg.o"
