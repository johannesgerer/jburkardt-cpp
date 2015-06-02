#!/bin/bash
#
cp filon.hpp /$HOME/include
#
g++ -c -I/$HOME/include filon.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling filon.cpp"
  exit
fi
rm compiler.txt
#
mv filon.o ~/libcpp/$ARCH/filon.o
#
echo "Library installed as ~/libcpp/$ARCH/filon.o"
