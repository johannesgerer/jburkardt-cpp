#!/bin/bash
#
cp square_arbq_rule.hpp /$HOME/include
#
g++ -c -I/$HOME/include square_arbq_rule.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling square_arbq_rule.cpp"
  exit
fi
#
mv square_arbq_rule.o ~/libcpp/$ARCH/square_arbq_rule.o
#
echo "Library installed as ~/libcpp/$ARCH/square_arbq_rule.o"
