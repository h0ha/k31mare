// u64mer.h -- packing char based 31-mers into uint64_t and back
//     2bit encoding used, only ACGT allowed, 64-bit flags error

#ifndef _U64MER_H_
#define _U64MER_H_

#include <cstdint>
#include <iostream>


// Converting from chars and packing into uint64_t

// bwa, Li's trick [ https://github.com/lh3/bwa ]
const uint64_t bwa_nt4_table[256] = {
     4, 4, 4, 4,              4, 4, 4, 4,     4, 4, 4, 4,     4, 4, 4, 4,
     4, 4, 4, 4,              4, 4, 4, 4,     4, 4, 4, 4,     4, 4, 4, 4,
     4, 4, 4, 4,              4, 4, 4, 4,     4, 4, 4, 4,     4, 4, 4, 4,
     4, 4, 4, 4,              4, 4, 4, 4,     4, 4, 4, 4,     4, 4, 4, 4,
     4, 0/*A*/, 4, 1/*C*/,    4, 4, 4, 2/*G*/,     4, 4, 4, 4,     4, 4, 4, 4,
     4, 4, 4, 4,              3/*T*/, 4, 4, 4,     4, 4, 4, 4,     4, 4, 4, 4,
     4, 0/*a*/, 4, 1/*c*/,    4, 4, 4, 2/*g*/,     4, 4, 4, 4,     4, 4, 4, 4,
     4, 4, 4, 4,              3/*t*/, 4, 4, 4,     4, 4, 4, 4,     4, 4, 4, 4,
     4, 4, 4, 4,              4, 4, 4, 4,     4, 4, 4, 4,     4, 4, 4, 4,
     4, 4, 4, 4,              4, 4, 4, 4,     4, 4, 4, 4,     4, 4, 4, 4,
     4, 4, 4, 4,              4, 4, 4, 4,     4, 4, 4, 4,     4, 4, 4, 4,
     4, 4, 4, 4,              4, 4, 4, 4,     4, 4, 4, 4,     4, 4, 4, 4,
     4, 4, 4, 4,              4, 4, 4, 4,     4, 4, 4, 4,     4, 4, 4, 4,
     4, 4, 4, 4,              4, 4, 4, 4,     4, 4, 4, 4,     4, 4, 4, 4,
     4, 4, 4, 4,              4, 4, 4, 4,     4, 4, 4, 4,     4, 4, 4, 4,
     4, 4, 4, 4,              4, 4, 4, 4,     4, 4, 4, 4,     4, 4, 4, 4
};



// packing using bit operators
static inline uint64_t u64mer (const char *kmer) {
    #define KM(pos) (bwa_nt4_table[kmer[pos]])
    #define KME(km) (KM(km) & 4)
    #define KMS(pos) (KM(pos) << (pos << 1))

    // python -c 'print(" | ".join(map(lambda x: "KME({})".format(x), range(31))))'
    // python -c 'print(" | ".join(map(lambda x: "KMS({})".format(x)), range(31))))'
    return ((
        KME(0) | KME(1) | KME(2) | KME(3) | KME(4) | KME(5) | KME(6) |
        KME(7) | KME(8) | KME(9) | KME(10) | KME(11) | KME(12) | KME(13) |
        KME(14) | KME(15) | KME(16) | KME(17) | KME(18) | KME(19) | KME(20) |
        KME(21) | KME(22) | KME(23) | KME(24) | KME(25) | KME(26) | KME(27) |
        KME(28) | KME(29) | KME(30)
        ) << 61 ) |
        ( 
            KMS(0) | KMS(1) | KMS(2) | KMS(3) | KMS(4) | KMS(5) | KMS(6) | KMS(7) | KMS(8) | KMS(9) | KMS(10) | KMS(11) | KMS(12) | KMS(13) | KMS(14) | KMS(15) | KMS(16) | KMS(17) | KMS(18) | KMS(19) | KMS(20) | KMS(21) | KMS(22) | KMS(23) | KMS(24) | KMS(25) | KMS(26) | KMS(27) | KMS(28) | KMS(29) | KMS(30)
        );

    #undef KME
    #undef KMS
    #undef KM
}

static inline uint64_t u64mer_2 (const char *kmer) {
    #define KM(pos) (bwa_nt4_table[kmer[pos]])
    #define KME(km) (KM(km) & 4)
    #define KMS(pos) (KM(pos) << (pos << 1))

    // python -c 'print(" | ".join(map(lambda x: "KME({})".format(x), range(31))))'
    // python -c 'print(" | ".join(map(lambda x: "KMS({})".format(x)), range(31))))'
    uint64_t err = (
        KME(0) | KME(1) | KME(2) | KME(3) | KME(4) | KME(5) | KME(6) |
        KME(7) | KME(8) | KME(9) | KME(10) | KME(11) | KME(12) | KME(13) |
        KME(14) | KME(15) | KME(16) | KME(17) | KME(18) | KME(19) | KME(20) |
        KME(21) | KME(22) | KME(23) | KME(24) | KME(25) | KME(26) | KME(27) |
        KME(28) | KME(29) | KME(30)
        );
    
    if (err)
        return (err << 61);

    return ( 
            KMS(0) | KMS(1) | KMS(2) | KMS(3) | KMS(4) | KMS(5) | KMS(6) | KMS(7) | KMS(8) | KMS(9) | KMS(10) | KMS(11) | KMS(12) | KMS(13) | KMS(14) | KMS(15) | KMS(16) | KMS(17) | KMS(18) | KMS(19) | KMS(20) | KMS(21) | KMS(22) | KMS(23) | KMS(24) | KMS(25) | KMS(26) | KMS(27) | KMS(28) | KMS(29) | KMS(30)
        );

    #undef KME
    #undef KMS
    #undef KM
}


// vectorized version of packing / never checked
static inline uint64_t u64mer_for (const char *kmer) {
    // res and err within same loop
    uint64_t res = 0;
    uint64_t err = 0;
    // may be use #pragma omp vector (static inline clash?) 
    for (auto i = 0; i < 31; i++) {
        uint64_t bits = bwa_nt4_table[kmer[i]];
        err |= bits & 4;
        res |= bits << (i << 1);
    }
    return (err << 61) | res;
}


// vectorized version of packing #2 / never checked
// separated error checking and kmer packing 
static inline uint64_t u64mer_for2 (const char *kmer) {
    // separate loop for res and err
    uint64_t res = 0;
    uint64_t err = 0;
    // may be use #pragma omp vector 
    for (auto i = 0; i < 31; i++) {
        err |= bwa_nt4_table[kmer[i]] & 4;
    }
    if (err)
        return (err << 61);
    // may be use #pragma omp vector 
    for (auto i = 0; i < 31; i++) {
        res |= bwa_nt4_table[kmer[i]] << (i << 1);
    }
    return res;
}

// test if uin64_t packed k-mer is a valid one
static inline bool test_u64mer(uint64_t val) {
   return !(val & (uint64_t(1) << 63));
}



// Unpacking from uint64_t into char[] buffer
const char u64mer_to_char[] = {'A', 'C', 'G', 'T'};

static inline void from_u64mer(uint64_t u64mer, char *buf, int k = 31) {
    if (buf == NULL)
        return;

    buf[0] = 0;
    buf[k] = 0;

    if (!test_u64mer) {
        memset(buf, 'N', k);
        return;
    }; 
    
    for (auto pos = 0; pos < k; ++pos) {
       uint8_t nt = (u64mer >> (pos << 1)) & 3;
       buf[pos] = u64mer_to_char[nt];
    }
    return;
}

#endif // _U64MER_H_
