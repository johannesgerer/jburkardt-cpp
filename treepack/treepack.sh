#!/bin/bash
#
cp treepack.hpp /$HOME/include
#
g++ -c -g -I/$HOME/include treepack.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling treepack.cpp"
  exit
fi
rm compiler.txt
#
mv treepack.o ~/libcpp/$ARCH/treepack.o
#
echo "Library installed as ~/libcpp/$ARCH/treepack.o"
