#!/bin/bash
#
cp asa183.hpp /$HOME/include
#
g++ -c -I /$HOME/include asa183.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling asa183.cpp"
  exit
fi
#
mv asa183.o ~/libcpp/$ARCH/asa183.o
#
echo "Library installed as ~/libcpp/$ARCH/asa183.o"
