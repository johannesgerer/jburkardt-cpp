#!/bin/bash
#
g++ -c -g -I/$HOME/include line_ncc_rule_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling line_ncc_rule_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ line_ncc_rule_prb.o /$HOME/libcpp/$ARCH/line_ncc_rule.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading line_ncc_rule_prb.o."
  exit
fi
#
rm line_ncc_rule_prb.o
#
mv a.out line_ncc_rule_prb
./line_ncc_rule_prb > line_ncc_rule_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running line_ncc_rule_prb."
  exit
fi
rm line_ncc_rule_prb
#
echo "Program output written to line_ncc_rule_prb_output.txt"
