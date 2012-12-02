#!/bin/bash
#
cp doomsday.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include doomsday.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling doomsday.cpp"
  exit
fi
rm compiler.txt
#
mv doomsday.o ~/libcpp/$ARCH/doomsday.o
#
echo "Library installed as ~/libcpp/$ARCH/doomsday.o"
