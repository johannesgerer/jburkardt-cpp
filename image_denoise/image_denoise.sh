#!/bin/bash
#
cp image_denoise.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include image_denoise.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling image_denoise.cpp"
  exit
fi
rm compiler.txt
#
mv image_denoise.o ~/libcpp/$ARCH/image_denoise.o
#
echo "Library installed as ~/libcpp/$ARCH/image_denoise.o"
