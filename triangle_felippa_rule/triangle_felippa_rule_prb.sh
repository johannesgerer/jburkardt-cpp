#!/bin/bash
#
g++ -c -I/$HOME/include triangle_felippa_rule_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle_felippa_rule_prb.cpp"
  exit
fi
#
g++ triangle_felippa_rule_prb.o /$HOME/libcpp/$ARCH/triangle_felippa_rule.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangle_felippa_rule_prb.o."
  exit
fi
#
rm triangle_felippa_rule_prb.o
#
mv a.out triangle_felippa_rule_prb
./triangle_felippa_rule_prb > triangle_felippa_rule_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running triangle_felippa_rule_prb."
  exit
fi
rm triangle_felippa_rule_prb
#
echo "Program output written to triangle_felippa_rule_prb_output.txt"
