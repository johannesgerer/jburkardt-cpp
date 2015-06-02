#!/bin/bash
#
g++ -c -I/$HOME/include tetrahedron_nco_rule_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling tetrahedron_nco_rule_prb.cpp"
  exit
fi
#
g++ tetrahedron_nco_rule_prb.o /$HOME/libcpp/$ARCH/tetrahedron_nco_rule.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading tetrahedron_nco_rule_prb.o."
  exit
fi
#
rm tetrahedron_nco_rule_prb.o
#
mv a.out tetrahedron_nco_rule_prb
./tetrahedron_nco_rule_prb > tetrahedron_nco_rule_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running tetrahedron_nco_rule_prb."
  exit
fi
rm tetrahedron_nco_rule_prb
#
echo "Program output written to tetrahedron_nco_rule_prb_output.txt"
