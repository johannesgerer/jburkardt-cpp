#!/bin/bash
#
cp ann.hpp ~/include
cp ann_perf.hpp ~/include
cp ann_x.hpp ~/include
cp bd_tree.hpp ~/include
cp kd_pr_search.hpp ~/include
cp kd_search.hpp ~/include
cp kd_tree.hpp ~/include
cp kd_util.hpp ~/include
cp pr_queue.hpp ~/include
cp pr_queue_k.hpp ~/include
#
g++ -c -O3 ann.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ann.cpp"
  exit
fi
rm compiler.txt
#
g++ -c -O3 bd_pr_search.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling bd_pr_search.cpp"
  exit
fi
rm compiler.txt
#
g++ -c -O3 bd_search.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling bd_search.cpp"
  exit
fi
rm compiler.txt
#
g++ -c -O3 bd_tree.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling bd_tree.cpp"
  exit
fi
rm compiler.txt
#
g++ -c -O3 brute.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling brute.cpp"
  exit
fi
rm compiler.txt
#
g++ -c -O3 kd_dump.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling kd_dump.cpp"
  exit
fi
rm compiler.txt
#
g++ -c -O3 kd_pr_search.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling kd_pr_search.cpp"
  exit
fi
rm compiler.txt
#
g++ -c -O3 kd_search.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling kd_search.cpp"
  exit
fi
rm compiler.txt
#
g++ -c -O3 kd_split.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling kd_split.cpp"
  exit
fi
rm compiler.txt
#
g++ -c -O3 kd_tree.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling kd_tree.cpp"
  exit
fi
rm compiler.txt
#
g++ -c -O3 kd_util.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling kd_util.cpp"
  exit
fi
rm compiler.txt
#
g++ -c -O3 perf.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling perf.cpp"
  exit
fi
rm compiler.txt
#
ar qc libann.a *.o
rm *.o
mv libann.a ~/libcpp/$ARCH/libann.a
#
echo "Library installed as ~/libcpp/$ARCH/libann.a"
