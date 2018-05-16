// processing.h -- main kmers processing functions utils

#ifndef _PROCESSING_H_
#define _PROCESSING_H_

#include <vector> 

#include "config.h"
#include "hash.h" // for fnv hashing
#include "io.h" // mmap and io stuff
#include "reads_kmers_stats.h" // reads, kmers stats
#include "ts.h" // timestamps

// using external library with fast and tightly packed hash map
// see https://github.com/sparsehash/sparsehash-c11
// using fnv hashing from GNU -- see hash.h
#include <sparsehash/sparse_hash_map>
typedef google::sparse_hash_map<uint64_t, counter_t, fnv_like_hash> seen_t;

// main file processing functions
// OpenMP threading and mmap-io is used in it
bool process_file(const FileInfo &file_info, const RunInfo &run_info,
         ReadsKmersStats &kmers_stats, TimeStamp &ts,
         std::vector<seen_t>& results
);

// fill k-mers into stat from the well formed chunk of reads
// assumed that start points to the proper read start and reads ids are starting from the same prefix of READID_PFX_LEN, see config.h
void fill_stats(char *start, char *end, seen_t &seen, ReadsKmersStats &kmers_stats, const seen_t* filter, char eol=EOL);

// util used to find the farest good read border
char *find_full_read_end(char *buf, size_t size, char eolc = EOL);


#endif // _PROCESSING_H_
