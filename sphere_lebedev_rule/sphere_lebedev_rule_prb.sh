#!/bin/bash
#
g++ -c -g -I/$HOME/include sphere_lebedev_rule_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sphere_lebedev_rule_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ sphere_lebedev_rule_prb.o /$HOME/libcpp/$ARCH/sphere_lebedev_rule.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sphere_lebedev_rule_prb.o."
  exit
fi
#
rm sphere_lebedev_rule_prb.o
#
mv a.out sphere_lebedev_rule_prb
./sphere_lebedev_rule_prb > sphere_lebedev_rule_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sphere_lebedev_rule_prb."
  exit
fi
rm sphere_lebedev_rule_prb
#
echo "Program output written to sphere_lebedev_rule_prb_output.txt"
