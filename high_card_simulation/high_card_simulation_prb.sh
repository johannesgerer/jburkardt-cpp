#!/bin/bash
#
g++ -c -I/$HOME/include high_card_simulation_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling high_card_simulation_prb.cpp"
  exit
fi
#
g++ high_card_simulation_prb.o /$HOME/libcpp/$ARCH/high_card_simulation.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading high_card_simulation_prb.o."
  exit
fi
#
rm high_card_simulation_prb.o
#
mv a.out high_card_simulation_prb
./high_card_simulation_prb > high_card_simulation_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running high_card_simulation_prb."
  exit
fi
rm high_card_simulation_prb
#
echo "Program output written to high_card_simulation_prb_output.txt"
