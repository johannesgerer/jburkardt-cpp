#!/bin/bash
#
cp change_making.hpp /$HOME/include
#
g++ -c -I/$HOME/include change_making.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling change_making.cpp"
  exit
fi
#
mv change_making.o ~/libcpp/$ARCH/change_making.o
#
echo "Library installed as ~/libcpp/$ARCH/change_making.o"
