#!/bin/bash
#
cp tetrahedron_arbq_rule.hpp /$HOME/include
#
g++ -c -I/$HOME/include tetrahedron_arbq_rule.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling tetrahedron_arbq_rule.cpp"
  exit
fi
#
mv tetrahedron_arbq_rule.o ~/libcpp/$ARCH/tetrahedron_arbq_rule.o
#
echo "Library installed as ~/libcpp/$ARCH/tetrahedron_arbq_rule.o"
