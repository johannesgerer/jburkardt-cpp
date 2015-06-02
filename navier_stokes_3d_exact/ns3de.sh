#!/bin/bash
#
cp ns3de.hpp /$HOME/include
#
g++ -c -I/$HOME/include ns3de.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling ns3de.cpp"
  exit
fi
#
mv ns3de.o ~/libcpp/$ARCH/ns3de.o
#
echo "Library installed as ~/libcpp/$ARCH/ns3de.o"
