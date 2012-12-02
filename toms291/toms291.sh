#!/bin/bash
#
cp toms291.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include toms291.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms291.cpp"
  exit
fi
rm compiler.txt
#
mv toms291.o ~/libcpp/$ARCH/toms291.o
#
echo "Library installed as ~/libcpp/$ARCH/toms291.o"
