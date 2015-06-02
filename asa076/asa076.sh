#!/bin/bash
#
cp asa076.hpp /$HOME/include
#
g++ -c -I /$HOME/include asa076.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling asa076.cpp"
  exit
fi
#
mv asa076.o ~/libcpp/$ARCH/asa076.o
#
echo "Library installed as ~/libcpp/$ARCH/asa076.o"
