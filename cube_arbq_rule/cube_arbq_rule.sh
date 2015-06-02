#!/bin/bash
#
cp cube_arbq_rule.hpp /$HOME/include
#
g++ -c -I/$HOME/include cube_arbq_rule.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling cube_arbq_rule.cpp"
  exit
fi
#
mv cube_arbq_rule.o ~/libcpp/$ARCH/cube_arbq_rule.o
#
echo "Library installed as ~/libcpp/$ARCH/cube_arbq_rule.o"
