#!/bin/bash
#
g++ -c -g truncated_normal_rule.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling truncated_normal_rule.cpp"
  exit
fi
#
g++ truncated_normal_rule.o -lm
if [ $? -ne 0 ]; then
  echo "Errors loading truncated_normal_rule.o"
  exit
fi
#
rm truncated_normal_rule.o
mv a.out ~/bincpp/$ARCH/truncated_normal_rule
#
echo "Executable installed as ~/bincpp/$ARCH/truncated_normal_rule."
