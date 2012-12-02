#!/bin/bash
#
g++ -c -g reactor_simulation.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling reactor_simulation.cpp"
  exit
fi
rm compiler.txt
#
g++ reactor_simulation.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading reactor_simulation.o"
  exit
fi
rm reactor_simulation.o
#
mv a.out reactor_simulation
./reactor_simulation > reactor_simulation_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running reactor_simulation"
  exit
fi
rm reactor_simulation
#
echo "Program output written to reactor_simulation_output.txt"
