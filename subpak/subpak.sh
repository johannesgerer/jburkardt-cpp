#!/bin/bash
#
cp subpak.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include subpak.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling subpak.cpp"
  exit
fi
rm compiler.txt
#
mv subpak.o ~/libcpp/$ARCH/subpak.o
#
echo "Library installed as ~/libcpp/$ARCH/subpak.o"
