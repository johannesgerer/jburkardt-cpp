#!/bin/bash
#
cp polpak.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include polpak.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling polpak.cpp"
  exit
fi
rm compiler.txt
#
mv polpak.o ~/libcpp/$ARCH/polpak.o
#
echo "Library installed as ~/libcpp/$ARCH/polpak.o"
