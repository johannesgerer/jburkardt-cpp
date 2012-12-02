#!/bin/bash
#
cp wandzura.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include wandzura.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling wandzura.cpp"
  exit
fi
rm compiler.txt
#
mv wandzura.o ~/libcpp/$ARCH/wandzura.o
#
echo "Library installed as ~/libcpp/$ARCH/wandzura.o"
