#!/bin/bash
#
cp kmeans.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include kmeans.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling kmeans.cpp"
  exit
fi
rm compiler.txt
#
mv kmeans.o ~/libcpp/$ARCH/kmeans.o
#
echo "Library installed as ~/libcpp/$ARCH/kmeans.o"
