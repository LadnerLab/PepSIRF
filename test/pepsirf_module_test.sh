#!/bin/bash

echo "PepSIRF Module Testing"
echo "======================"

echo "[Normalize Module]"

let passed=0
let failed=0
#
echo "Combination: diff approach, a single matrix, name filter"
echo "-------"
../build/pepsirf norm \
-p 'test_raw_matrix.tsv' \
-a 'diff' \
-n 'test_samplenames_norm.txt' \
--precision '8' \
-o 'norm_module_testing_02.10/test_diff_names_norm_output.tsv'
echo "-------"
#
echo "Combination: ratio approach, a single matrix, name filter"
echo "-------"
../build/pepsirf norm \
-p 'test_raw_matrix.tsv' \
-a 'ratio' \
-n 'test_samplenames_norm.txt' \
--precision '8' \
-o 'norm_module_testing_02.10/test_ratio_names_norm_output.tsv'
echo "-------"
#
echo "Combination: diff ratio approach, a single matrix, name filter"
echo "-------"
../build/pepsirf norm \
-p 'test_raw_matrix.tsv' \
-a 'diff_ratio' \
-n 'test_samplenames_norm.txt' \
--precision '8' \
-o 'norm_module_testing_02.10/test_diffratio_names_norm_output.tsv'
echo "-------"
#
echo "Combination: diff approach, two matrices, name filter"
echo "-------"
../build/pepsirf norm \
-p 'test_raw_matrix.tsv' \
--negative_control 'test_raw_matrix.tsv' \
-n 'test_samplenames_norm.txt' \
-a 'diff' \
--precision '8' \
-o 'norm_module_testing_02.10/test_diff_separate_neg_control_names_norm_output.tsv'
echo "-------"
#
echo "Combination: ratio approach, two matrices, name filter"
echo "-------"
../build/pepsirf norm \
-p 'test_raw_matrix.tsv' \
--negative_control 'test_raw_matrix.tsv' \
-n 'test_samplenames_norm.txt' \
-a 'ratio' \
--precision '8' \
-o 'norm_module_testing_02.10/test_ratio_separate_neg_control_names_norm_output.tsv'
echo "-------"
#
echo "Combination: diff ratio approach, two matrices, name filter"
echo "-------"
../build/pepsirf norm \
-p 'test_raw_matrix.tsv' \
--negative_control 'test_raw_matrix.tsv' \
-n 'test_samplenames_norm.txt' \
-a 'diff_ratio' \
--precision '8' \
-o 'norm_module_testing_02.10/test_diffratio_separate_neg_control_names_norm_output.tsv'
echo "-------"
#
echo "Combination: diff approach, single matrix, id filter"
echo "-------"
../build/pepsirf norm \
-p 'test_raw_matrix.tsv' \
-a 'diff' \
-s 'S_' \
--precision '8' \
-o 'norm_module_testing_02.10/test_diff_id_norm_output.tsv'
echo "-------"
#
echo "Combination: ratio approach, single matrix, id filter"
echo "-------"
../build/pepsirf norm \
-p 'test_raw_matrix.tsv' \
-a 'ratio' \
-s 'S_' \
--precision '8' \
-o 'norm_module_testing_02.10/test_ratio_id_norm_output.tsv'
echo "-------"
#
echo "Combination: diff ratio approach, single matrix, id filter"
echo "-------"
../build/pepsirf norm \
-p 'test_raw_matrix.tsv' \
-a 'diff_ratio' \
-s 'S_' \
--precision '8' \
-o 'norm_module_testing_02.10/test_diffratio_id_norm_output.tsv'
echo "-------"
#
echo "Combination: diff approach, two matrices, id filter"
echo "-------"
../build/pepsirf norm \
-p 'test_raw_matrix.tsv' \
-a 'diff' \
--negative_control 'test_raw_matrix.tsv' \
-s 'S_' \
--precision '8' \
-o 'norm_module_testing_02.10/test_diff_separate_neg_control_id_norm_output.tsv'
echo "-------"
#
echo "Combination: ratio approach, two matrices, id filter"
echo "-------"
../build/pepsirf norm \
-p 'test_raw_matrix.tsv' \
-a 'ratio' \
--negative_control 'test_raw_matrix.tsv' \
-s 'S_' \
--precision '8' \
-o 'norm_module_testing_02.10/test_ratio_separate_neg_control_id_norm_output.tsv'
echo "-------"
#
echo "Combination: diff ratio approach, two matrices, id filter"
echo "-------"
../build/pepsirf norm \
-p 'test_raw_matrix.tsv' \
-a 'diff_ratio' \
--negative_control 'test_raw_matrix.tsv' \
-s 'S_' \
--precision '8' \
-o 'norm_module_testing_02.10/test_diffratio_separate_neg_control_id_norm_output.tsv'
echo "-------"
#
echo "Combination: size factors approach, single matrix"
echo "-------"
../build/pepsirf norm \
-p 'test_raw_matrix.tsv' \
-a 'size_factors' \
--precision '8' \
-o 'norm_module_testing_02.10/test_diff_names_norm_output.tsv'
echo "-------"
#
echo "Combination: col sum approach, single matrix"
echo "-------"
../build/pepsirf norm \
-p 'test_raw_matrix.tsv' \
-a 'col_sum' \
--precision '8' \
-o 'norm_module_testing_02.10/test_diff_names_norm_output.tsv'
echo "-------"


