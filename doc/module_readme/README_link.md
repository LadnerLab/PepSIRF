#### Module Link
```
PepSIRF: Peptide-based Serological Immune Response Framework link module.
:
  -h [ --help ]                         Produce help message and exit.
                                        
  -o [ --output ] arg (=link_output.tsv)
                                        Name of the file to write output to. 
                                        Output will be in the form of a 
                                        tab-delimited file with a header. Each 
                                        entry will be of the form:
                                        peptide_name TAB id:score,id:score, and
                                        so on. 
                                        
  --protein_file arg                    Name of fasta file containing protein 
                                        sequences from which a design was 
                                        created.
                                        
  --peptide_file arg                    Name of fasta file containing aa 
                                        peptides that have been designed as 
                                        part of a library.
                                        
  --tax_id_index arg (=1)               The index (0-based, valid values 
                                        include 0-3) of the tax id to use for 
                                        linking peptides and species. For 
                                        example, if this argument is passed 
                                        with the value 1, 
                                        the species ID will be used. (2 for 
                                        genus, 3 for family. 0 can vary 
                                        depending upon the 
                                        method used for assigning the 0'th ID.
                                        
  --kmer_redundancy_control             Control for kmer redundancy when 
                                        creating the peptide linkage map. 
                                        Instead of a peptide receiving one 
                                        point for each kmer it receives for a 
                                        species, it recieves 1 / ( the number 
                                        of times the kmer appears in the 
                                        original design ) points.
                                        
  -k [ --kmer_size ] arg                Kmer size to use when creating the 
                                        linkage map. Entries in the linkage 
                                        file will contain peptides and the 
                                        species ids of species that share a 
                                        kmer with this peptide. For example, if
                                        k is 7 and there exists a line in the 
                                        linkage file of the form: 
                                         'peptide_1 TAB 455:12,423:10'
                                         then peptide_1 shares 12 7-mers with 
                                        the species with id '455', and 10 
                                        7-mers with the species that has id 
                                        423.
```
