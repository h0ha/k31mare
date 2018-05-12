// ts.h -- timestamps and stuff

#ifndef _TS_H
#define _TS_H

#include <omp.h>
#include <iostream>

struct TimeStamp {
    double epoch = get_wtime();
    double epoch_end = 0.0;

    double in_iterations = 0.0;
    double in_iterations_all = 0.0;
    double in_mmap_io = 0.0;
    double in_processing = 0.0;
    double in_merge = 0.0;
    double in_filter = 0.0;

    double _start_io = 0;
    double _start_iteration = 0;
    double _start_iterations_all = 0;
    double _start_processing = 0;
    double _start_merge = 0;
    double _start_filter = 0;


    inline double get_wtime() {
        return omp_get_wtime();
    }

    // global start / end
    inline void start_epoch() {
        epoch = get_wtime();    
    }
    inline void end_epoch() {
        epoch_end = get_wtime();    
    }

    // whole processing loop
    inline void start_iterations() {
        _start_iterations_all = get_wtime();
    } 

    inline void end_iterations() {
       in_iterations_all += get_wtime() - _start_iterations_all;
    }

    // iterations inside loop
    inline void start_iteration() {
        _start_iteration = get_wtime();    
    }

    inline void end_iteration() {
       in_iterations += get_wtime() - _start_iteration;
    }

    // mmap io related stuff
    inline void start_io() {
        _start_io = get_wtime();    
    }

    inline void end_io() {
       in_mmap_io += get_wtime() - _start_io;
    }

    // kmers processing
    inline void start_processing() {
        _start_processing = get_wtime();    
    }

    inline void end_processing() {
       in_processing += get_wtime() - _start_processing;
    }


    // merge
    inline void start_merge() {
        _start_merge = get_wtime();
    } 

    inline void end_merge() {
       in_merge += get_wtime() - _start_merge;
    }

    // filter
    inline void start_filter() {
        _start_filter = get_wtime();
    } 

    inline void end_filter() {
       in_filter += get_wtime() - _start_filter;
    }

    friend std::ostream& operator<<(std::ostream& out, const TimeStamp& ts) {
        out << "#T total: " << ts.epoch_end - ts.epoch << std::endl;
        out << "#T alloc: " << ts._start_iterations_all - ts.epoch << std::endl;
        out << "#T all iterations: " << ts.in_iterations_all << std::endl;
        out << "#T in iterations(master): " << ts.in_iterations << std::endl;
        out << "#T processing(master): " << ts.in_processing << std::endl;
        out << "#T map_io(master): " << ts.in_mmap_io << std::endl;
        out << "#T hash merging: " << ts.in_merge << std::endl;
        out << "#T hash filtering: " << ts.in_filter << std::endl;

        return out;
    }
};

#endif // _TS_H
