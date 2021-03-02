#!/bin/bash

# How to use:
# Please provide the following space delimited arguments in this order:
# $1= build directory
# $2= output directory
# $3= threshfile filename
# $4= sample pairs list
# $5= raw scores matrix filename
# $6= raw score constraint

echo "PepSIRF P_enrich Module Testing"
echo "=================================================================="

mkdir $2
# 16 directories created, each will contain output for a given combo
mkdir $2'/1_matrix_single_thresh_no_raw'
mkdir $2'/1_matrix_dual_thresh_no_raw'
mkdir $2'/2_matrix_single_thresh_no_raw'
mkdir $2'/2_matrix_dual_thresh_no_raw'
mkdir $2'/2_matrix_combo_thresh_no_raw'
mkdir $2'/3_matrix_single_thresh_no_raw'
mkdir $2'/3_matrix_dual_thresh_no_raw'
mkdir $2'/3_matrix_combo_thresh_no_raw'
mkdir $2'/1_matrix_single_thresh_raw'
mkdir $2'/1_matrix_dual_thresh_raw'
mkdir $2'/2_matrix_single_thresh_raw'
mkdir $2'/2_matrix_dual_thresh_raw'
mkdir $2'/2_matrix_combo_thresh_raw'
mkdir $2'/3_matrix_single_thresh_raw'
mkdir $2'/3_matrix_dual_thresh_raw'
mkdir $2'/3_matrix_combo_thresh_raw'


# Create a single matrix, single threshold file
data="$(echo "$a" | sed '1!d' $3)" # gives first line
IFS=$'\t' read -ra matrix <<< $data # get matrix filename
IFS=$',' read -ra thresh <<< ${matrix[1]} # get thresholds
echo ${matrix[0]}$'\t'${thresh[1]} >> threshfile_1matrix_1thresh.txt
echo "Thresholds input file 'threshfile_1matrix_1thresh.txt' created."
echo ${matrix[0]}$'\t'${thresh[0]}','${thresh[1]} >> threshfile_1matrix_2thresh.txt
echo "Thresholds input file 'threshfile_1matrix_2thresh.txt' created."

echo ${matrix[0]}$'\t'${thresh[1]} >> threshfile_2matrix_1thresh.txt
echo "Thresholds input file 'threshfile_2matrix_1thresh.txt' created."
echo ${matrix[0]}$'\t'${thresh[1]} >> threshfile_3matrix_1thresh.txt
echo "Thresholds input file 'threshfile_3matrix_1thresh.txt' created."
echo ${matrix[0]}$'\t'${thresh[0]}','${thresh[1]} >> threshfile_2matrix_2thresh.txt
echo "Thresholds input file 'threshfile_2matrix_2thresh.txt' created."
echo ${matrix[0]}$'\t'${thresh[0]}','${thresh[1]} >> threshfile_3matrix_2thresh.txt
echo "Thresholds input file 'threshfile_3matrix_2thresh.txt' created."
echo ${matrix[0]}$'\t'${thresh[1]} >> threshfile_2matrix_combo.txt
echo "Thresholds input file 'threshfile_2matrix_comb.txt' created."
echo ${matrix[0]}$'\t'${thresh[0]}','${thresh[1]} >> threshfile_3matrix_combo.txt
echo "Thresholds input file 'threshfile_3matrix_comb.txt' created."

data="$(echo "$a" | sed '2!d' $3)"
IFS=$'\t' read -ra matrix <<< $data
IFS=$',' read -ra thresh <<< ${matrix[1]}
echo ${matrix[0]}$'\t'${thresh[1]} >> threshfile_2matrix_1thresh.txt
echo ${matrix[0]}$'\t'${thresh[0]}','${thresh[1]} >> threshfile_2matrix_2thresh.txt
echo ${matrix[0]}$'\t'${thresh[0]}','${thresh[1]} >> threshfile_2matrix_combo.txt

echo ${matrix[0]}$'\t'${thresh[1]} >> threshfile_3matrix_1thresh.txt
echo ${matrix[0]}$'\t'${thresh[0]}','${thresh[1]} >> threshfile_3matrix_2thresh.txt
echo ${matrix[0]}$'\t'${thresh[1]} >> threshfile_3matrix_combo.txt

data="$(echo "$a" | sed '3!d' $3)"
IFS=$'\t' read -ra matrix <<< $data
IFS=$',' read -ra thresh <<< ${matrix[1]}
echo ${matrix[0]}$'\t'${thresh[1]} >> threshfile_3matrix_1thresh.txt
echo ${matrix[0]}$'\t'${thresh[0]}','${thresh[1]} >> threshfile_3matrix_2thresh.txt
echo ${matrix[0]}$'\t'${thresh[0]}','${thresh[1]} >> threshfile_3matrix_combo.txt

echo "Threshold files have successfully been created."

echo "Combination: 1 matrix file, single threshold, no raw score"
echo "----------------"
$1/pepsirf p_enrich \
-t threshfile_1matrix_1thresh.txt \
-s $4 \
-o $2'/1_matrix_single_thresh_no_raw'
echo '----------------------------------------------------------'

echo "Combination: 1 matrix file, dual thresholds, no raw score"
echo "----------------"
$1/pepsirf p_enrich \
-t threshfile_1matrix_2thresh.txt \
-s $4 \
-o $2'/1_matrix_dual_thresh_no_raw'
echo '----------------------------------------------------------'

echo "Combination: 2 matrix file, single threshold, no raw score"
echo "----------------"
$1/pepsirf p_enrich \
-t threshfile_2matrix_1thresh.txt \
-s $4 \
-o $2'/2_matrix_single_thresh_no_raw'
echo '----------------------------------------------------------'

echo "Combination: 2 matrix file, dual thresholds, no raw score"
echo "----------------"
$1/pepsirf p_enrich \
-t threshfile_2matrix_2thresh.txt \
-s $4 \
-o $2'/2_matrix_dual_thresh_no_raw'
echo '----------------------------------------------------------'

echo "Combination: 2 matrix file, threshold combo, no raw score"
echo "----------------"
$1/pepsirf p_enrich \
-t threshfile_2matrix_combo.txt \
-s $4 \
-o $2'/2_matrix_combo_thresh_no_raw'
echo '----------------------------------------------------------'

echo "Combination: 3 matrix file, single threshold, no raw score"
echo "----------------"
$1/pepsirf p_enrich \
-t threshfile_3matrix_1thresh.txt \
-s $4 \
-o $2'/3_matrix_single_thresh_no_raw'
echo '----------------------------------------------------------'

echo "Combination: 3 matrix file, dual thresholds, no raw score"
echo "----------------"
$1/pepsirf p_enrich \
-t threshfile_3matrix_2thresh.txt \
-s $4 \
-o $2'/3_matrix_dual_thresh_no_raw'
echo '----------------------------------------------------------'

echo "Combination: 3 matrix file, threshold combo, no raw score"
echo "----------------"
$1/pepsirf p_enrich \
-t threshfile_3matrix_combo.txt \
-s $4 \
-o $2'/3_matrix_combo_thresh_no_raw'
echo '----------------------------------------------------------'

echo "Combination: 1 matrix file, single threshold, raw score"
echo "----------------"
$1/pepsirf p_enrich \
-t threshfile_1matrix_1thresh.txt \
-s $4 \
-r $5 \
--raw_score_constraint $6 \
-o $2'/1_matrix_single_thresh_raw'
echo '----------------------------------------------------------'

echo "Combination: 1 matrix file, dual thresholds, raw score"
echo "----------------"
$1/pepsirf p_enrich \
-t threshfile_1matrix_2thresh.txt \
-s $4 \
-r $5 \
--raw_score_constraint $6 \
-o $2'/1_matrix_dual_thresh_raw'
echo '----------------------------------------------------------'

echo "Combination: 2 matrix file, single threshold, raw score"
echo "----------------"
$1/pepsirf p_enrich \
-t threshfile_2matrix_1thresh.txt \
-s $4 \
-r $5 \
--raw_score_constraint $6 \
-o $2'/2_matrix_single_thresh_raw'
echo '----------------------------------------------------------'

echo "Combination: 2 matrix file, dual thresholds, raw score"
echo "----------------"
$1/pepsirf p_enrich \
-t threshfile_2matrix_2thresh.txt \
-s $4 \
-r $5 \
--raw_score_constraint $6 \
-o $2'/2_matrix_dual_thresh_raw'
echo '----------------------------------------------------------'

echo "Combination: 2 matrix file, threshold combo, raw score"
echo "----------------"
$1/pepsirf p_enrich \
-t threshfile_2matrix_combo.txt \
-s $4 \
-r $5 \
--raw_score_constraint $6 \
-o $2'/2_matrix_combo_thresh_raw'
echo '----------------------------------------------------------'

echo "Combination: 3 matrix file, single threshold, raw score"
echo "----------------"
$1/pepsirf p_enrich \
-t threshfile_3matrix_1thresh.txt \
-s $4 \
-r $5 \
--raw_score_constraint $6 \
-o $2'/3_matrix_single_thresh_raw'
echo '----------------------------------------------------------'

echo "Combination: 3 matrix file, dual thresholds, raw score"
echo "----------------"
$1/pepsirf p_enrich \
-t threshfile_3matrix_2thresh.txt \
-s $4 \
-r $5 \
--raw_score_constraint $6 \
-o $2'/3_matrix_dual_thresh_raw'
echo '----------------------------------------------------------'

echo "Combination: 3 matrix file, threshold combo, raw score"
echo "----------------"
$1/pepsirf p_enrich \
-t threshfile_3matrix_combo.txt \
-s $4 \
-r $5 \
--raw_score_constraint $6 \
-o $2'/3_matrix_combo_thresh_raw'
echo '----------------------------------------------------------'

echo "End of P_enrich Testing"
echo "=================================================================="

echo "Comparing P_enrich output to Existing True Output"
echo "----------------------------------------------------------"
diff -yr $2'/1_matrix_single_thresh_no_raw' 'true_penrich_output/1_matrix_single_thresh_no_raw'
diff -yr $2'/1_matrix_dual_thresh_no_raw' 'true_penrich_output/1_matrix_dual_thresh_no_raw'
diff -yr $2'/2_matrix_single_thresh_no_raw' 'true_penrich_output/2_matrix_single_thresh_no_raw'
diff -yr $2'/2_matrix_dual_thresh_no_raw' 'true_penrich_output/2_matrix_dual_thresh_no_raw'
diff -yr $2'/2_matrix_combo_thresh_no_raw' 'true_penrich_output/2_matrix_combo_thresh_no_raw'
diff -yr $2'/3_matrix_single_thresh_no_raw' 'true_penrich_output/3_matrix_single_thresh_no_raw'
diff -yr $2'/3_matrix_dual_thresh_no_raw' 'true_penrich_output/3_matrix_dual_thresh_no_raw'
diff -yr $2'/3_matrix_combo_thresh_no_raw' 'true_penrich_output/3_matrix_combo_thresh_no_raw'
diff -yr $2'/1_matrix_single_thresh_raw' 'true_penrich_output/1_matrix_single_thresh_raw'
diff -yr $2'/1_matrix_dual_thresh_raw' 'true_penrich_output/1_matrix_dual_thresh_raw'
diff -yr $2'/2_matrix_single_thresh_raw' 'true_penrich_output/2_matrix_single_thresh_raw'
diff -yr $2'/2_matrix_dual_thresh_raw' 'true_penrich_output/2_matrix_dual_thresh_raw'
diff -yr $2'/2_matrix_combo_thresh_raw' 'true_penrich_output/2_matrix_combo_thresh_raw'
diff -yr $2'/3_matrix_single_thresh_raw' 'true_penrich_output/3_matrix_single_thresh_raw'
diff -yr $2'/3_matrix_dual_thresh_raw' 'true_penrich_output/3_matrix_dual_thresh_raw'
diff -yr $2'/3_matrix_combo_thresh_raw' 'true_penrich_output/3_matrix_combo_thresh_raw'

diff -sr $2'/1_matrix_single_thresh_no_raw' 'true_penrich_output/1_matrix_single_thresh_no_raw'
diff -sr $2'/1_matrix_dual_thresh_no_raw' 'true_penrich_output/1_matrix_dual_thresh_no_raw'
diff -sr $2'/2_matrix_single_thresh_no_raw' 'true_penrich_output/2_matrix_single_thresh_no_raw'
diff -sr $2'/2_matrix_dual_thresh_no_raw' 'true_penrich_output/2_matrix_dual_thresh_no_raw'
diff -sr $2'/2_matrix_combo_thresh_no_raw' 'true_penrich_output/2_matrix_combo_thresh_no_raw'
diff -sr $2'/3_matrix_single_thresh_no_raw' 'true_penrich_output/3_matrix_single_thresh_no_raw'
diff -sr $2'/3_matrix_dual_thresh_no_raw' 'true_penrich_output/3_matrix_dual_thresh_no_raw'
diff -sr $2'/3_matrix_combo_thresh_no_raw' 'true_penrich_output/3_matrix_combo_thresh_no_raw'
diff -sr $2'/1_matrix_single_thresh_raw' 'true_penrich_output/1_matrix_single_thresh_raw'
diff -sr $2'/1_matrix_dual_thresh_raw' 'true_penrich_output/1_matrix_dual_thresh_raw'
diff -sr $2'/2_matrix_single_thresh_raw' 'true_penrich_output/2_matrix_single_thresh_raw'
diff -sr $2'/2_matrix_dual_thresh_raw' 'true_penrich_output/2_matrix_dual_thresh_raw'
diff -sr $2'/2_matrix_combo_thresh_raw' 'true_penrich_output/2_matrix_combo_thresh_raw'
diff -sr $2'/3_matrix_single_thresh_raw' 'true_penrich_output/3_matrix_single_thresh_raw'
diff -sr $2'/3_matrix_dual_thresh_raw' 'true_penrich_output/3_matrix_dual_thresh_raw'
diff -sr $2'/3_matrix_combo_thresh_raw' 'true_penrich_output/3_matrix_combo_thresh_raw'
echo "----------------------------------------------------------"
echo " "
echo "End of P_enrich Testing"
echo "=================================================================="
