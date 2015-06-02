#!/bin/bash
#
g++ -c -I/$HOME/include line_felippa_rule_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling line_felippa_rule_prb.cpp"
  exit
fi
#
g++ line_felippa_rule_prb.o /$HOME/libcpp/$ARCH/line_felippa_rule.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading line_felippa_rule_prb.o."
  exit
fi
#
rm line_felippa_rule_prb.o
#
mv a.out line_felippa_rule_prb
./line_felippa_rule_prb > line_felippa_rule_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running line_felippa_rule_prb."
  exit
fi
rm line_felippa_rule_prb
#
echo "Program output written to line_felippa_rule_prb_output.txt"
