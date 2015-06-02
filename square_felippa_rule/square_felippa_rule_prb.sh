#!/bin/bash
#
g++ -c -I/$HOME/include square_felippa_rule_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling square_felippa_rule_prb.cpp"
  exit
fi
#
g++ square_felippa_rule_prb.o /$HOME/libcpp/$ARCH/square_felippa_rule.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading square_felippa_rule_prb.o."
  exit
fi
#
rm square_felippa_rule_prb.o
#
mv a.out square_felippa_rule_prb
./square_felippa_rule_prb > square_felippa_rule_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running square_felippa_rule_prb."
  exit
fi
rm square_felippa_rule_prb
#
echo "Program output written to square_felippa_rule_prb_output.txt"
