#include "peptide_bin.h"




// bin_collection 
bin_collection::iterator bin_collection::begin() { return bins.begin(); }
bin_collection::iterator bin_collection::end() { return bins.end(); }

bin_collection::const_iterator bin_collection::begin() const { return bins.begin(); }
bin_collection::const_iterator bin_collection::end() const { return bins.end(); }


