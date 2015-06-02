#!/bin/bash
#
g++ htmlindex.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling htmlindex.cpp"
  exit
fi
#
mv a.out ~/bincpp/$ARCH/htmlindex
#
echo "Program installed as ~/bincpp/$ARCH/htmlindex."
