// hash.h -- wrapping GNU FNV hashing
//     see FNV-hashing http://isthe.com/chongo/tech/comp/fnv/#FNV-param

#ifndef _HASH_H
#define _HASH_H

#include <cstdint>
#include <bits/hash_bytes.h>


struct fnv_like_hash
    {
        inline size_t operator() (uint64_t __x) const
        {
            return std::_Fnv_hash_bytes((void*)&__x, 8, 16777619UL);
        }
    };


// as proposed by the author, one can use this approach for reduced bit-sizes
struct fnv_like_hash_28
    {
        inline size_t operator() (uint64_t __x) const
        {
            size_t hash = std::_Fnv_hash_bytes((void*)&__x, 8, 16777619UL);
            return ((hash>>28) ^ hash) & ((1UL << 28) - 1);
        }
    };


#endif // _HASH_H
