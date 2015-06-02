#!/bin/bash
#
cp asa066.hpp /$HOME/include
#
g++ -c -I /$HOME/include asa066.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling asa066.cpp"
  exit
fi
#
mv asa066.o ~/libcpp/$ARCH/asa066.o
#
echo "Library installed as ~/libcpp/$ARCH/asa066.o"
