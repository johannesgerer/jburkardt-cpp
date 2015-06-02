#!/bin/bash
#
cp tetrahedron_keast_rule.hpp /$HOME/include
#
g++ -c -I /$HOME/include tetrahedron_keast_rule.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling tetrahedron_keast_rule.cpp"
  exit
fi
#
mv tetrahedron_keast_rule.o ~/libcpp/$ARCH/tetrahedron_keast_rule.o
#
echo "Library installed as ~/libcpp/$ARCH/tetrahedron_keast_rule.o"
