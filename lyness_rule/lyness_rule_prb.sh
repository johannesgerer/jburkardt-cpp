#!/bin/bash
#
g++ -c -g -I/$HOME/include lyness_rule_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling lyness_rule_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ lyness_rule_prb.o /$HOME/libcpp/$ARCH/lyness_rule.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading lyness_rule_prb.o."
  exit
fi
#
rm lyness_rule_prb.o
#
mv a.out lyness_rule_prb
./lyness_rule_prb > lyness_rule_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running lyness_rule_prb."
  exit
fi
rm lyness_rule_prb
#
echo "Program output written to lyness_rule_prb_output.txt"
