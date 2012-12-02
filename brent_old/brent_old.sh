#!/bin/bash
#
cp brent_old.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include brent_old.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling brent_old.cpp"
  exit
fi
rm compiler.txt
#
mv brent_old.o ~/libcpp/$ARCH/brent_old.o
#
echo "Library installed as ~/libcpp/$ARCH/brent_old.o"
