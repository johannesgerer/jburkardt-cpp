#!/bin/bash
#
cp chrpak.hpp /$HOME/include
#
g++ -c -I /$HOME/include chrpak.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling chrpak.cpp"
  exit
fi
#
mv chrpak.o ~/libcpp/$ARCH/chrpak.o
#
echo "Library installed as ~/libcpp/$ARCH/chrpak.o"
