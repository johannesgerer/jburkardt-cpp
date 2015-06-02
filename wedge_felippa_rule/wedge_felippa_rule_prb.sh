#!/bin/bash
#
g++ -c -I/$HOME/include wedge_felippa_rule_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling wedge_felippa_rule_prb.cpp"
  exit
fi
#
g++ -o wedge_felippa_rule_prb wedge_felippa_rule_prb.o /$HOME/libcpp/$ARCH/wedge_felippa_rule.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading wedge_felippa_rule_prb.o."
  exit
fi
#
rm wedge_felippa_rule_prb.o
#
./wedge_felippa_rule_prb > wedge_felippa_rule_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running wedge_felippa_rule_prb."
  exit
fi
rm wedge_felippa_rule_prb
#
echo "Program output written to wedge_felippa_rule_prb_output.txt"
