// reads_kmers_stats.h -- reads and kmers statistics storage
#ifndef  _READS_KMERS_STATS_H
#define  _READS_KMERS_STATS_H

#include <iostream>

struct ReadsKmersStats {
    size_t total_reads = 0;
    size_t total_kmers = 0;
    size_t bad_kmers = 0;


    void clear() {
        total_reads = 0;
        total_kmers = 0;
        bad_kmers = 0;
    }

    friend std::ostream& operator<<(std::ostream& out, const ReadsKmersStats& s) {
        out << "total reads: " << s.total_reads << std::endl << std::flush;
        out << "bad kmers: " << s.bad_kmers
            << "  of " << s.total_kmers
            << " ( " << size_t((1000.0L * s.bad_kmers) / s.total_kmers) / 10 << "% ) "
            << std::endl << std::flush;
        return out;
    }

    ReadsKmersStats& operator += (const ReadsKmersStats &other) {
        total_reads += other.total_reads;
        total_kmers += other.total_kmers;
        bad_kmers += other.bad_kmers;
        return *this;
    }
};

#endif // _READS_KMERS_STATS_H
