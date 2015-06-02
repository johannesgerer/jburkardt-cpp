#!/bin/bash
#
g++ -c -I/$HOME/include cube_felippa_rule_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling cube_felippa_rule_prb.cpp"
  exit
fi
#
g++ cube_felippa_rule_prb.o /$HOME/libcpp/$ARCH/cube_felippa_rule.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cube_felippa_rule_prb.o."
  exit
fi
#
rm cube_felippa_rule_prb.o
#
mv a.out cube_felippa_rule_prb
./cube_felippa_rule_prb > cube_felippa_rule_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running cube_felippa_rule_prb."
  exit
fi
rm cube_felippa_rule_prb
#
echo "Program output written to cube_felippa_rule_prb_output.txt"
