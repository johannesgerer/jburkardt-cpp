#!/bin/bash
#
cp pink_noise.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include pink_noise.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pink_noise.cpp"
  exit
fi
rm compiler.txt
#
mv pink_noise.o ~/libcpp/$ARCH/pink_noise.o
#
echo "Library installed as ~/libcpp/$ARCH/pink_noise.o"
