#!/bin/bash
#
cp tetrahedron_ncc_rule.hpp /$HOME/include
#
g++ -c -I /$HOME/include tetrahedron_ncc_rule.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling tetrahedron_ncc_rule.cpp"
  exit
fi
#
mv tetrahedron_ncc_rule.o ~/libcpp/$ARCH/tetrahedron_ncc_rule.o
#
echo "Library installed as ~/libcpp/$ARCH/tetrahedron_ncc_rule.o"
