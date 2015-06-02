#!/bin/bash
#
cp asa144.hpp /$HOME/include
#
g++ -c -I /$HOME/include asa144.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling asa144.cpp"
  exit
fi
#
mv asa144.o ~/libcpp/$ARCH/asa144.o
#
echo "Library installed as ~/libcpp/$ARCH/asa144.o"
