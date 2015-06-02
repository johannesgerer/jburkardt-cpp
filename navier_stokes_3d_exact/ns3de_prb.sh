#!/bin/bash
#
g++ -c -I/$HOME/include ns3de_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling ns3de_prb.cpp"
  exit
fi
#
g++ -o ns3de_prb ns3de_prb.o /$HOME/libcpp/$ARCH/ns3de.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ns3de_prb.o."
  exit
fi
#
rm ns3de_prb.o
#
./ns3de_prb > ns3de_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ns3de_prb."
  exit
fi
rm ns3de_prb
#
echo "Program output written to ns3de_prb_output.txt"
