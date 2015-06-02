#!/bin/bash
#
g++ -c -g -I/$HOME/include circle_rule_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling circle_rule_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ circle_rule_prb.o /$HOME/libcpp/$ARCH/circle_rule.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading circle_rule_prb.o."
  exit
fi
#
rm circle_rule_prb.o
#
mv a.out circle_rule_prb
./circle_rule_prb > circle_rule_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running circle_rule_prb."
  exit
fi
rm circle_rule_prb
#
echo "Program output written to circle_rule_prb_output.txt"
