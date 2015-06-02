#!/bin/bash
#
cp asa007.hpp /$HOME/include
#
g++ -c -I /$HOME/include asa007.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling asa007.cpp"
  exit
fi
#
mv asa007.o ~/libcpp/$ARCH/asa007.o
#
echo "Library installed as ~/libcpp/$ARCH/asa007.o"
