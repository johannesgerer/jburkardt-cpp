#!/bin/bash
#
cp naca.hpp /$HOME/include
#
g++ -c -I/$HOME/include naca.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling naca.cpp"
  exit
fi
rm compiler.txt
#
mv naca.o ~/libcpp/$ARCH/naca.o
#
echo "Library installed as ~/libcpp/$ARCH/naca.o"
