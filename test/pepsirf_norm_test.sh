#!/bin/bash

# How to use:
# Please provide the following space delimited arguments in this order:
# $1= build directory
# $2= output directory
# $3= raw score matrix filename
# $4= list of names for filtering in norm
# $5= id for filtering in norm
# example
# ./pepsirf_norm_test.sh ~/Documents/work/PepSIRF/build ~/Documents/work/PepSIRF/test/norm_penrich_module_testing ~/Documents/work/PepSIRF/test/test_raw_matrix.tsv S_000,S_001,S_002,S_003,S_004,S_005,S_006,S_007 S_

echo "PepSIRF Norm Module Testing"
echo "=================================================================="

mkdir $2
mkdir $2'/norm_single_matrix_output'
mkdir $2'/norm_dual_matrix_output'
mkdir $2'/norm_python_output'

echo "Combination: diff approach, a single matrix, name filter"
echo "----------------"
$1/pepsirf norm \
-p $3 \
--precision 4 \
-a 'diff' \
-n $4 \
-o $2'/norm_single_matrix_output/test_diff_names_norm_output.tsv'
echo '----------------------------------------------------------'



echo "Combination: ratio approach, a single matrix, name filter"
echo "----------------"
$1/pepsirf norm \
-p $3 \
--precision 4 \
-a 'ratio' \
-n $4 \
-o $2'/norm_single_matrix_output/test_ratio_names_norm_output.tsv'
echo "----------------------------------------------------------"



echo "Combination: diff ratio approach, a single matrix, name filter"
echo "----------------"
$1/pepsirf norm \
-p $3 \
--precision 4 \
-a 'diff_ratio' \
-n $4 \
-o $2'/norm_single_matrix_output/test_diffratio_names_norm_output.tsv'
echo "----------------------------------------------------------"



echo "Combination: diff approach, two matrices, name filter"
echo "----------------"
$1/pepsirf norm \
-p $3 \
--precision 4 \
--negative_control $3 \
-n $4 \
-a 'diff' \
-o $2'/norm_dual_matrix_output/test_diff_names_norm_output.tsv'
echo "----------------------------------------------------------"



echo "Combination: ratio approach, two matrices, name filter"
echo "----------------"
$1/pepsirf norm \
-p $3 \
--precision 4 \
--negative_control $3 \
-n $4 \
-a 'ratio' \
-o $2'/norm_dual_matrix_output/test_ratio_names_norm_output.tsv'
echo "----------------------------------------------------------"



echo "Combination: diff ratio approach, two matrices, name filter"
echo "----------------"
$1/pepsirf norm \
-p $3 \
--precision 4 \
--negative_control $3 \
-n $4 \
-a 'diff_ratio' \
-o $2'/norm_dual_matrix_output/test_diffratio_names_norm_output.tsv'
echo "----------------------------------------------------------"



echo "Combination: diff approach, single matrix, id filter"
echo "----------------"
$1/pepsirf norm \
-p $3 \
--precision 4 \
-a 'diff' \
-s $5 \
-o $2'/norm_single_matrix_output/test_diff_id_norm_output.tsv'
echo "----------------------------------------------------------"



echo "Combination: ratio approach, single matrix, id filter"
echo "----------------"
$1/pepsirf norm \
-p $3 \
--precision 4 \
-a 'ratio' \
-s $5 \
-o $2'/norm_single_matrix_output/test_ratio_id_norm_output.tsv'
echo "----------------------------------------------------------"



echo "Combination: diff ratio approach, single matrix, id filter"
echo "----------------"
$1/pepsirf norm \
-p $3 \
--precision 4 \
-a 'diff_ratio' \
-s $5 \
-o $2'/norm_single_matrix_output/test_diffratio_id_norm_output.tsv'
echo "----------------------------------------------------------"



echo "Combination: diff approach, two matrices, id filter"
echo "----------------"
$1/pepsirf norm \
-p $3 \
--precision 4 \
-a 'diff' \
--negative_control $3 \
-s $5 \
-o $2'/norm_dual_matrix_output/test_diff_id_norm_output.tsv'
echo "----------------------------------------------------------"



echo "Combination: ratio approach, two matrices, id filter"
echo "----------------"
$1/pepsirf norm \
-p $3 \
--precision 4 \
-a 'ratio' \
--negative_control $3 \
-s $5 \
-o $2'/norm_dual_matrix_output/test_ratio_id_norm_output.tsv'
echo "----------------------------------------------------------"



echo "Combination: diff ratio approach, two matrices, id filter"
echo "----------------"
$1/pepsirf norm \
-p $3 \
--precision 4 \
-a 'diff_ratio' \
--negative_control $3 \
-s $5 \
-o $2'/norm_dual_matrix_output/test_diffratio_id_norm_output.tsv'
echo "----------------------------------------------------------"



echo "Combination: size factors approach, single matrix"
echo "----------------"
$1/pepsirf norm \
-p $3 \
--precision 4 \
-a 'size_factors' \
-o $2'/test_size_factors_norm_output.tsv'
echo "----------------------------------------------------------"



echo "Combination: col sum approach, single matrix"
echo "----------------"
$1/pepsirf norm \
-p $3 \
--precision 4 \
-a 'col_sum' \
-o $2'/test_col_sum_norm_output.tsv'
echo "----------------------------------------------------------"
echo " "
echo "Comparing Single Matrix to Dual Matrix"
echo "----------------------------------------------------------"
diff -yr $2'/norm_dual_matrix_output' $2'/norm_single_matrix_output' > $2'/norm_matrix_usage_comparison_results.txt'
diff -sr $2'/norm_dual_matrix_output' $2'/norm_single_matrix_output'
echo "----------------------------------------------------------"

echo " "
echo "Comparing ID to Namelist"
echo "----------------------------------------------------------"
diff -s $2'/norm_single_matrix_output/test_diff_names_norm_output.tsv' $2'/norm_single_matrix_output/test_diff_id_norm_output.tsv'
diff -s $2'/norm_single_matrix_output/test_ratio_names_norm_output.tsv' $2'/norm_single_matrix_output/test_ratio_id_norm_output.tsv'
diff -s $2'/norm_single_matrix_output/test_diffratio_names_norm_output.tsv' $2'/norm_single_matrix_output/test_diffratio_id_norm_output.tsv'
diff -s $2'/norm_dual_matrix_output/test_diff_names_norm_output.tsv' $2'/norm_dual_matrix_output/test_diff_id_norm_output.tsv'
diff -s $2'/norm_dual_matrix_output/test_ratio_names_norm_output.tsv' $2'/norm_dual_matrix_output/test_ratio_id_norm_output.tsv'
diff -s $2'/norm_dual_matrix_output/test_diffratio_names_norm_output.tsv' $2'/norm_dual_matrix_output/test_diffratio_id_norm_output.tsv'
echo "----------------------------------------------------------"
echo " "

python3 $1/../extensions/SBnorm.py \
-i $3 \
-a 'diff' \
-n $4 \
-o $2'/norm_python_output/test_diff_names_norm_output.tsv'

python3 $1/../extensions/SBnorm.py \
-i $3 \
-a 'ratio' \
-n $4 \
-o $2'/norm_python_output/test_ratio_names_norm_output.tsv'

python3 $1/../extensions/SBnorm.py \
-i $3 \
-a 'diffratio' \
-n $4 \
-o $2'/norm_python_output/test_diffratio_names_norm_output.tsv'

echo "Testing PepSIRF Results Against Python Results"

python3 $1/../extensions/compMatrices.py \
-1 $2'/norm_python_output/test_diff_names_norm_output.tsv' \
-2 $2'/norm_dual_matrix_output/test_diff_names_norm_output.tsv'

python3 $1/../extensions/compMatrices.py \
-1 $2'/norm_python_output/test_ratio_names_norm_output.tsv' \
-2 $2'/norm_dual_matrix_output/test_ratio_names_norm_output.tsv'

python3 $1/../extensions/compMatrices.py \
-1 $2'/norm_python_output/test_diffratio_names_norm_output.tsv' \
-2 $2'/norm_dual_matrix_output/test_diffratio_names_norm_output.tsv'
echo " "
echo "End of Norm Testing"
echo "=================================================================="
