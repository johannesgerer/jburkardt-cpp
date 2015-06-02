#!/bin/bash
#
cp square_symq_rule.hpp /$HOME/include
#
g++ -c -I/$HOME/include square_symq_rule.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling square_symq_rule.cpp"
  exit
fi
#
mv square_symq_rule.o ~/libcpp/$ARCH/square_symq_rule.o
#
echo "Library installed as ~/libcpp/$ARCH/square_symq_rule.o"
