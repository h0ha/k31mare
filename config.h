// config.h -- various constants and types

#ifndef _CONFIG_H
#define _CONFIG_H

#include <cstdint>
#include <vector>

// just for info, k=31 hardcoded, see u64mer.h
#define KMER_LEN (31)

// kmer-counts default storage size. set to 120*1e6 if have lot's of mem
#define DEFAULT_DICT_SIZE (40*1e6)


// kmer-counts type and tuple type used for final median sort and output steps
typedef uint32_t counter_t;
typedef std::tuple<uint64_t, counter_t, std::vector<counter_t>> tuple_t;


// unix EOL
#define EOL '\n'


// mmap sizes and mlock usage
#define BUFFER_SIZE_IN_PAGES (256)

// mlocking is safer for multithreading,
// but be sure to have enough ulimits for this
#define USE_MLOCK true


// do not run in parallel if have fewer data
#define MIN_DATA_SIZE_IN_PAGES (2)

// if fill_stats should check for the unexpected read end/start
#define CHECK_READ_START
#define READID_PFX_LEN (4)

// print out info to std::cerr each number of steps
#define PRINT_INFO_EACH (500)

#endif //_CONFIG_H
