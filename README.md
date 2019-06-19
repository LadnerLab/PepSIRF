# PepSIRF: Peptide-based Serological Immune Response Framework  


### Dependencies
- cmake version 3.9 or higher.
- boost version 3.6 or higher.

### Compiling/Building 
```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make pep_sirf
```

Alternatively, one can use the included build script.
```
chmod +x build.sh
./build.sh
```
Both of these options will create a ```pep_sirf``` executable in the build directory.

### Running Tests
```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make pepsirf_test

./pepsirf_test
```
#### Note: Building for HPC Clusters
Some additional arguments may need to be specified, depending upon how your system's
module system is managed. For example on NAU's Monsoon cluster:
```
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DBOOST_ROOT=/packages/boost/1.66.0-gcc-6.2.0 ..
```

### Usage
PepSIRF is a module-based program, and each module has its own set of arguments. 
```
USAGE: pep_sirf [ --help | module_name <module_args*> ] 
--help, -h displays this message, while 'pep_sirf module_name --help' will display the help for the module module_name.
```

Currently only forward indexes and runs where samples have already been indexed on the Reverse indexes are supported.

#### Module Demux
```
PepSIRF: Peptide-based Serological Immune Response Framework demultiplexing module. 
This module takes the following parameters and outputs the counts of each reference 
sequence for each sample.:
  -h [ --help ]                        Produce help message
  --input_r1 arg                       Input forward reads fastq file to parse.
                                       
  --index arg                          Name of fasta file containing forward 
                                       and (if included )reverse index 
                                       sequences.
                                       
  -l [ --library ] arg                 Designed library containing nucleic acid
                                       peptides. Library should be in fasta 
                                       format and should contain sequences that
                                       were used to design input_r1.
                                       
  -r [ --read_per_loop ] arg (=100000) The number of fastq records read a time.
                                       A higher value will result in more 
                                       memory usage by the program, but will 
                                       also result in fewer disk accesses, 
                                       increasing performance of the program.
                                       
  -m [ --max_mismatches ] arg (=0)     The maximum number of 'mismatches' to 
                                       tolerate when parsing reads. If a read 
                                       is not within this value of any of the 
                                       sequences within the designed library it
                                       will not be considered. Note that here 
                                       we define a mismatch by the Hamming 
                                       distance D between a reference sequence 
                                       r and a read sequence s. If D( r, s ) <=
                                       max_mismatches we say that r and s are 
                                       similar. Note that if for some read q D(
                                       r, s ) <= max_mismatches and D( r, q ) 
                                       <= max_mismatches, then we say that the 
                                       sequence whose distance is the minimum 
                                       between D( r, s ) and D( r, q ) maps to 
                                       reference r. 
                                       
  --seq_start arg                      Start index (0-based) of each read where
                                       we expect the designed peptide to begin.
                                       For each read, we start at this index 
                                       and read for seq_len number of 
                                       characters. Remember: this argument 
                                       should be zero-based!
                                       
  --seq_len arg                        The length of the designed peptides. 
                                       Note that we assume all of the designed 
                                       peptides are the same length.
                                       
  --f_index_start arg                  Start index (0-based) of each read where
                                       we expect the forward index sequences to
                                       be found.
                                       
  --f_index_len arg                    Length of forward index sequences. For 
                                       each read we start at f_index_start and 
                                       grab f_index_len nucleotides.
                                       
  -o [ --output ] arg (=output.csv)    The name of the output file to write 
                                       counts to. Each line in this file will 
                                       be a tab-separated list of values, 
                                       where each entry i is either the name of
                                       a sequence or the counts for this 
                                       sequence in sample i. This file will 
                                       have a header labelling each column, 
                                       i'th tab-separated value of column i 
                                       of the header will be the sample name of
                                       sample i. If we traverse this column, we
                                       will see the count of this sample for 
                                       each sequence. 
                                       
  -s [ --samplelist ] arg              A tab-delimited list of samples, one 
                                       sample per line. If the samples are 
                                       already indexed by I2 only the forward 
                                       index (I1) and the sample name are 
                                       required. The first item in each 
                                       tab-delimited line is the forward (I1) 
                                       index, the second (if included) is the 
                                       reverse (I2) index, and the third is the
                                       samplename. 
                                       
  -t [ --num_threads ] arg (=2)        Number of threads to use for analyses.
```