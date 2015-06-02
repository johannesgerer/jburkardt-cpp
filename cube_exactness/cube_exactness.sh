#!/bin/bash
#
cp cube_exactness.hpp /$HOME/include
#
g++ -c -I/$HOME/include cube_exactness.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling cube_exactness.cpp"
  exit
fi
#
mv cube_exactness.o ~/libcpp/$ARCH/cube_exactness.o
#
echo "Library installed as ~/libcpp/$ARCH/cube_exactness.o"
