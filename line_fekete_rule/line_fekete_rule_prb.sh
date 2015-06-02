#!/bin/bash
#
g++ -c -I/$HOME/include line_fekete_rule_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling line_fekete_rule_prb.cpp"
  exit
fi
#
g++ line_fekete_rule_prb.o /$HOME/libcpp/$ARCH/line_fekete_rule.o /$HOME/libcpp/$ARCH/qr_solve.o /$HOME/libcpp/$ARCH/r8lib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading line_fekete_rule_prb.o."
  exit
fi
#
rm line_fekete_rule_prb.o
#
mv a.out line_fekete_rule_prb
./line_fekete_rule_prb > line_fekete_rule_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running line_fekete_rule_prb."
  exit
fi
rm line_fekete_rule_prb
#
echo "Program output written to line_fekete_rule_prb_output.txt"
