#!/bin/bash
#
cp c4lib.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include c4lib.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling c4lib.cpp"
  exit
fi
rm compiler.txt
#
mv c4lib.o ~/libcpp/$ARCH/c4lib.o
#
echo "Library installed as ~/libcpp/$ARCH/c4lib.o"
