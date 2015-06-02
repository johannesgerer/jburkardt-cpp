#!/bin/bash
#
cp cube_felippa_rule.hpp /$HOME/include
#
g++ -c -I /$HOME/include cube_felippa_rule.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling cube_felippa_rule.cpp"
  exit
fi
#
mv cube_felippa_rule.o ~/libcpp/$ARCH/cube_felippa_rule.o
#
echo "Library installed as ~/libcpp/$ARCH/cube_felippa_rule.o"
