#!/bin/bash
#
cp image_edge.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include image_edge.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling image_edge.cpp"
  exit
fi
rm compiler.txt
#
mv image_edge.o ~/libcpp/$ARCH/image_edge.o
#
echo "Library installed as ~/libcpp/$ARCH/image_edge.o"
