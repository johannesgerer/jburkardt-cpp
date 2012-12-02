#!/bin/bash
#
cp compass_search.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include compass_search.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling compass_search.cpp."
  exit
fi
rm compiler.txt
#
mv compass_search.o ~/libcpp/$ARCH/compass_search.o
#
echo "Library installed as ~/libcpp/$ARCH/compass_search.o"
