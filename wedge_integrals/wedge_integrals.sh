#!/bin/bash
#
cp wedge_integrals.hpp /$HOME/include
#
g++ -c -I/$HOME/include wedge_integrals.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling wedge_integrals.cpp"
  exit
fi
#
mv wedge_integrals.o ~/libcpp/$ARCH/wedge_integrals.o
#
echo "Library installed as ~/libcpp/$ARCH/wedge_integrals.o"
