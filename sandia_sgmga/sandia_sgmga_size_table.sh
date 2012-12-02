#!/bin/bash
#
g++ -c -g -I/$HOME/include sandia_sgmga_size_table.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sandia_sgmga_size_table.cpp"
  exit
fi
rm compiler.txt
#
g++ sandia_sgmga_size_table.o /$HOME/libcpp/$ARCH/sandia_sgmga.o /$HOME/libcpp/$ARCH/sandia_rules2.o /$HOME/libcpp/$ARCH/sandia_rules.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sandia_sgmga_size_table.o."
  exit
fi
#
rm sandia_sgmga_size_table.o
#
mv a.out sandia_sgmga_size_table
./sandia_sgmga_size_table > sandia_sgmga_size_table_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sandia_sgmga_size_table."
  exit
fi
rm sandia_sgmga_size_table
#
echo "Program output written to sandia_sgmga_size_table_output.txt"
