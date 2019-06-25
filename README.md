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
sequence for each sample. Later we reference the 'distance' between sequences. For this module we define distance by the Hamming distance D between a reference sequence r and a read sequence s. If D( r, s ) <= max_mismatches we say that r and s are similar. Note that if for some read q D( r, s ) <= max_mismatches and D( r, q ) <= max_mismatches, then we say that the sequence whose distance is the minimum between D( r, s ) and D( r, q ) maps to reference r. Additionally if D( r, s ) == D( r, q ) then we discard the read as we cannot say whether r maps to s or q. 
:
  -h [ --help ]                        Produce help message
  --input_r1 arg                       Input forward reads fastq file to parse.
                                       
  --input_r2 arg                       Input reverse reads fastq file to parse.
                                       Note that if this argument is not 
                                       supplied only forward indices will be 
                                       used to identify samples.
                                       
  --index arg                          Name of fasta file containing forward 
                                       and (potentially) reverse index 
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
                                       
  --f_index arg                        Positional values for f_index. This 
                                       argument must be passed as 3 
                                       comma-separated values. The first item 
                                       represents the (0-based) expected start 
                                       index of the forward index. The second 
                                       represents the length of the forward 
                                       index, and the third represents the 
                                       number of mismatches that are tolerated 
                                       for this index. An example is '--f_index
                                       12,12,2'. This says that we start at 
                                       (0-based) index 12, grab the next 12 
                                       characters, and if a perfect match is 
                                       not found for these grabbed characters 
                                       we look for a match to the forward index
                                       sequences with up to two allowed 
                                       mismatches.
                                       .
  --seq arg                            Positional values for nucleotide 
                                       sequence data. This argument must be 
                                       passed as 3 comma-separated values. The 
                                       first item represents the (0-based) 
                                       expected start index of the sequence. 
                                       The second represents the length of the 
                                       sequence, and the third represents the 
                                       number of mismatches that are tolerated 
                                       for a sequence. An example is '--seq 
                                       43,90,2'. This says that we start at 
                                       (0-based) index 43, grab the next 90 
                                       characters, and if a perfect match is 
                                       not found for these grabbed characters 
                                       we look for a match to the designed 
                                       library sequences with up to two allowed
                                       mismatches.
                                       .
  --r_index arg (=0,0,0)               Positional values for r_index. This 
                                       argument must be passed as 3 
                                       comma-separated values. The first item 
                                       represents the (0-based) expected start 
                                       index of the reverse index. The second 
                                       represents the length of the reverse 
                                       index, and the third represents the 
                                       number of mismatches that are tolerated 
                                       for this index. An example is '--r_index
                                       12,12,2'. This says that we start at 
                                       (0-based) index 12, grab the next 12 
                                       characters, and if a perfect match is 
                                       not found for these grabbed characters 
                                       we look for a match to the reverse index
                                       sequences with up to two allowed 
                                       mismatches.
                                       .
  --concatemer arg                     Concatenated primer sequences. If this 
                                       concatemer is found within a read, we 
                                       know that a potential sequence from the 
                                       designed library was not included. The 
                                       number of times this concatemer is 
                                       recorded in the input file is reported.
                                       
  -o [ --output ] arg (=output.csv)    The name of the output file to write 
                                       counts to. Each line in this file will 
                                       be a tab-separated list of values, where
                                       each entry i is either the name of a 
                                       sequence or the counts for this sequence
                                       in sample i. This file will have a 
                                       header labelling each column, i'th 
                                       tab-separated value of column i of the 
                                       header will be the sample name of sample
                                       i. If we traverse this column, we will 
                                       see the count of this sample for each 
                                       sequence. 
                                       
  -a [ --aa_counts ] arg               The name of the file to write aggregated
                                       counts to when aa sequences from a 
                                       designed library have multiple different
                                       nt encodings. If this option is 
                                       included, names of sequences in the file
                                       supplied by the --library flag MUST be 
                                       of the form ID-NUM, where ID can contain
                                       any characters but '-', and NUM 
                                       represents the id of this encoding. ID 
                                       and NUM MUST be separated by a single 
                                       dash '-' character. If supplied 
                                       aggregated counts for each sequence will
                                       be written to this file. For example, 
                                       suppose we have TG1_1-1 and TG1_1-2 in 
                                       our library, which says that we 
                                       generated two encodings for TG1_1. If 
                                       included, the file will have a single 
                                       TG1_1 entry, where the count in column i
                                       is the sum of the value of column i from
                                       TG1_1-1 and TG1_1-2.
                                       
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

### Generating Code Documentation
To create documentation for the code run the following:
```
cd doc
doxygen doxygen.conf
```
From here you can either use the html or LaTeX versions.
If you choose to use the LaTeX version, run the following:
```
cd latex
make
```
This creates a file named ```refman.pdf```.
