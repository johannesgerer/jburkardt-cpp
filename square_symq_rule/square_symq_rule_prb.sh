#!/bin/bash
#
g++ -c -I/$HOME/include square_symq_rule_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling square_symq_rule_prb.cpp"
  exit
fi
#
g++ -o square_symq_rule_prb square_symq_rule_prb.o /$HOME/libcpp/$ARCH/square_symq_rule.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading square_symq_rule_prb.o."
  exit
fi
#
rm square_symq_rule_prb.o
#
./square_symq_rule_prb > square_symq_rule_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running square_symq_rule_prb."
  exit
fi
rm square_symq_rule_prb
#
echo "Program output written to square_symq_rule_prb_output.txt"
