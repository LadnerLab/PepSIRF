#### Module Demux
```
PepSIRF: Peptide-based Serological Immune Response Framework demultiplexing module. 
This module takes the following parameters and outputs the counts of each reference 
sequence for each sample. Later we reference the 'distance' between sequences. For this module we define distance by the Hamming distance D between a reference sequence r and a read sequence s. If D( r, s ) <= max_mismatches we say that r and s are similar. Note that if for some read q D( r, s ) <= max_mismatches and D( r, q ) <= max_mismatches, then we say that the sequence whose distance is the minimum between D( r, s ) and D( r, q ) maps to reference r. Additionally if D( r, s ) == D( r, q ) then we discard the read as we cannot say whether r maps to s or q. 
:
  -h [ --help ]                        Produce help message
  --input_r1 arg                       Input forward reads fastq file to parse.
                                       
  --input_r2 arg                       Input reverse reads fastq file to parse. Note that if this argument is not supplied only forward indices will be used to identify samples.
                                       
  -i [ --index ] arg                   Name of fasta file containing forward and (potentially) reverse index sequences.
                                       
  -l [ --library ] arg                 Designed library containing nucleic acid peptides. Library should be in fasta format and should contain sequences that were used to design input_r1 If
                                       this flag is not included, reference-independent demultiplexing will be performed. In this mode, each probe in the region specified by '--seq' will be
                                       considered its own reference. The count of each probe will be the number of times it has appeared.
                                       
  -r [ --read_per_loop ] arg (=100000) The number of fastq records read a time. A higher value will result in more memory usage by the program, but will also result in fewer disk accesses, 
                                       increasing performance of the program.
                                       
  --f_index arg                        Positional values for f_index. This argument must be passed as 3 comma-separated values. The first item represents the (0-based) expected start index 
                                       of the forward index. The second represents the length of the forward index, and the third represents the number of mismatches that are tolerated for 
                                       this index. An example is '--f_index 12,12,2'. This says that we start at (0-based) index 12, grab the next 12 characters, and if a perfect match is 
                                       not found for these grabbed characters we look for a match to the forward index sequences with up to two allowed mismatches.
                                       .
  --seq arg                            Positional values for nucleotide sequence data. This argument must be passed as 3 comma-separated values. The first item represents the (0-based) 
                                       expected start index of the sequence. The second represents the length of the sequence, and the third represents the number of mismatches that are 
                                       tolerated for a sequence. An example is '--seq 43,90,2'. This says that we start at (0-based) index 43, grab the next 90 characters, and if a perfect 
                                       match is not found for these grabbed characters we look for a match to the designed library sequences with up to two allowed mismatches.
                                       .
  --r_index arg (=0,0,0)               Positional values for r_index. This argument must be passed as 3 comma-separated values. The first item represents the (0-based) expected start index 
                                       of the reverse index. The second represents the length of the reverse index, and the third represents the number of mismatches that are tolerated for 
                                       this index. An example is '--r_index 12,12,2'. This says that we start at (0-based) index 12, grab the next 12 characters, and if a perfect match is 
                                       not found for these grabbed characters we look for a match to the reverse index sequences with up to two allowed mismatches.
                                       .
  -c [ --concatemer ] arg              Concatenated primer sequences. If this concatemer is found within a read, we know that a potential sequence from the designed library was not 
                                       included. The number of times this concatemer is recorded in the input file is reported.
                                       
  -o [ --output ] arg (=output.csv)    The name of the output file to write counts to. Each line in this file will be a tab-separated list of values, where each entry i is either the name 
                                       of a sequence or the counts for this sequence in sample i. This file will have a header labelling each column, i'th tab-separated value of column i of
                                       the header will be the sample name of sample i. If we traverse this column, we will see the count of this sample for each sequence. 
                                       
  -a [ --aa_counts ] arg               The name of the file to write aggregated counts to when aa sequences from a designed library have multiple different nt encodings. If this option is 
                                       included, names of sequences in the file supplied by the --library flag MUST be of the form ID-NUM, where ID can contain any characters but '-', and 
                                       NUM represents the id of this encoding. ID and NUM MUST be separated by a single dash '-' character. If supplied aggregated counts for each sequence 
                                       will be written to this file. For example, suppose we have TG1_1-1 and TG1_1-2 in our library, which says that we generated two encodings for TG1_1. 
                                       If included, the file will have a single TG1_1 entry, where the count in column i is the sum of the value of column i from TG1_1-1 and TG1_1-2.
                                       
  -s [ --samplelist ] arg              A tab-delimited list of samples, one sample per line. If the samples are already indexed by I2 only the forward index (I1) and the sample name are 
                                       required. The first item in each tab-delimited line is the forward (I1) index, the second (if included) is the reverse (I2) index, and the third is 
                                       the samplename. 
                                       
  --phred_base arg (=33)               Phred base to use when parsing fastq quality scores. Valid options include 33 or 64.
                                       
  --phred_min_score arg (=0)           The minimum average score for the sequence portion of a read for it to be considered for matching. This means that if the average phred33/64 score for
                                       a read at the expected locations of a library sequence is not at least this value then the read will be discarded.
                                       
  -t [ --num_threads ] arg (=2)        Number of threads to use for analyses.
```