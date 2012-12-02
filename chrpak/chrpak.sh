#!/bin/bash
#
cp chrpak.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include chrpak.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling chrpak.cpp"
  exit
fi
rm compiler.txt
#
mv chrpak.o ~/libcpp/$ARCH/chrpak.o
#
echo "Library installed as ~/libcpp/$ARCH/chrpak.o"
