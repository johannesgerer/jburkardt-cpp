#!/bin/bash
#
#  Make fundamental library.
#
echo "Making FSU_SUB."
bash fsu_sub.bash
if [ $? -ne 0 ]; then
  echo "Errors in fsu_sub.bash."
  exit
fi
#
#  Make other librarys.
#
echo "Making FSU_CVT, _HALTON, _HAMMERSLEY, _LATINIZE, _QUALITY."
bash fsu_cvt.bash
if [ $? -ne 0 ]; then
  echo "Errors in fsu_cvt.bash."
  exit
fi
bash fsu_halton.bash
if [ $? -ne 0 ]; then
  echo "Errors in fsu_halton.bash."
  exit
fi
bash fsu_hammersley.bash
if [ $? -ne 0 ]; then
  echo "Errors in fsu_hammersley.bash."
  exit
fi
bash fsu_latinize.bash
if [ $? -ne 0 ]; then
  echo "Errors in fsu_latinize.bash."
  exit
fi
bash fsu_quality.bash
if [ $? -ne 0 ]; then
  echo "Errors in fsu_quality.bash."
  exit
fi
#
#  Make standalone codes
#
echo "Making FSU_CVT_STANDALONE, _HALTON_, _HAMMERSLEY_, _LATINIZE_, _QUALITY_"
bash fsu_cvt_standalone.bash
if [ $? -ne 0 ]; then
  echo "Errors in fsu_cvt_standalone.bash."
  exit
fi
bash fsu_halton_standalone.bash
if [ $? -ne 0 ]; then
  echo "Errors in fsu_halton_standalone.bash."
  exit
fi
bash fsu_hammersley_standalone.bash
if [ $? -ne 0 ]; then
  echo "Errors in fsu_hammersley_standalone.bash."
  exit
fi
bash fsu_latinize_standalone.bash
if [ $? -ne 0 ]; then
  echo "Errors in fsu_latinize_standalone.bash."
  exit
fi
bash fsu_quality_standalone.bash
if [ $? -ne 0 ]; then
  echo "Errors in fsu_quality_standalone.bash."
  exit
fi
#
#  Run problem sets.
#
echo "Making FSU_CVT_PROBLEMS, _LATINIZE_, _QUALITY_."
bash fsu_cvt_problems.bash
if [ $? -ne 0 ]; then
  echo "Errors in fsu_cvt_problems.bash."
  exit
fi
bash fsu_latinize_problems.bash
if [ $? -ne 0 ]; then
  echo "Errors in fsu_latinize_problems.bash."
  exit
fi
bash fsu_quality_problems.bash
if [ $? -ne 0 ]; then
  echo "Errors in fsu_quality_problems.bash."
  exit
fi
#
echo "Normal end of execution."

