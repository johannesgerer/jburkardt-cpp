#!/bin/bash
#
g++ -c -g -I/$HOME/include gm_rule_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling gm_rule_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ gm_rule_prb.o /$HOME/libcpp/$ARCH/gm_rule.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading gm_rule_prb.o."
  exit
fi
#
rm gm_rule_prb.o
#
mv a.out gm_rule_prb
./gm_rule_prb > gm_rule_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running gm_rule_prb."
  exit
fi
rm gm_rule_prb
#
echo "Program output written to gm_rule_prb_output.txt"
