#!/bin/bash
#
cp polpak.hpp /$HOME/include
#
g++ -c -I /$HOME/include polpak.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling polpak.cpp"
  exit
fi
#
mv polpak.o ~/libcpp/$ARCH/polpak.o
#
echo "Library installed as ~/libcpp/$ARCH/polpak.o"
