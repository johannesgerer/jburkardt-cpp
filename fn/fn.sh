#!/bin/bash
#
cp fn.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include fn.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fn.cpp"
  exit
fi
rm compiler.txt
#
mv fn.o ~/libcpp/$ARCH/fn.o
#
echo "Library installed as ~/libcpp/$ARCH/fn.o"
