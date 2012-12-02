#!/bin/bash
#
g++ -c -g -I/$HOME/include sgmg_size_table.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sgmg_size_table.cpp"
  exit
fi
rm compiler.txt
#
g++ sgmg_size_table.o /$HOME/libcpp/$ARCH/sgmg.o /$HOME/libcpp/$ARCH/sandia_rules.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sgmg_size_table.o."
  exit
fi
#
rm sgmg_size_table.o
#
mv a.out sgmg_size_table
./sgmg_size_table > sgmg_size_table_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sgmg_size_table."
  exit
fi
rm sgmg_size_table
#
echo "Program output written to sgmg_size_table_output.txt"
