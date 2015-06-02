#!/bin/bash
#
cp simplex_gm_rule.hpp /$HOME/include
#
g++ -c -I /$HOME/include simplex_gm_rule.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling simplex_gm_rule.cpp"
  exit
fi
#
mv simplex_gm_rule.o ~/libcpp/$ARCH/simplex_gm_rule.o
#
echo "Library installed as ~/libcpp/$ARCH/simplex_gm_rule.o"
