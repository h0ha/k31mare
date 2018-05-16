// k31mare.cpp -- main process

#include <iostream>
#include <queue>
#include <tuple>

#include <omp.h>

#include "config.h"
#include "io.h"
#include "processing.h"
#include "reads_kmers_stats.h"
#include "u64mer.h"


// print usage
void print_usage_err(char *argv[]) {
    std::cerr << "usage: " << argv[0] << " [--no-size-sort] in_file1.fastq [ in_file2.fastq ... ]" << std::endl << std::flush;
}

// main
int main(int argc, char* argv[]) {
    // failing if not enough arguments
    if (argc < 2) {
        std::cerr << "no input files specified" <<  std::endl << std::flush;
        print_usage_err(argv);
        return 1;
    }

    // checking if parameter disabling sorting by filesize is set
    int argv_files_start = 1;
    if (std::string(argv[1]) == std::string("--no-size-sort"))
        ++argv_files_start;

    // preparin file info structures with names, sizes, fds
    std::vector<FileInfo> files;
    for (auto i = argv_files_start; i < argc; i++) {
        FileInfo fi(argv[i]);
        if (fi) files.push_back(std::move(fi));
    }

    // process files, smaller first -- hope to filter not-final k-mers
    if (std::string(argv[1]) != std::string("--no-size-sort")) {
        std::sort(files.begin(), files.end(), [](const auto& l, const auto& r) {
            return l.file_size < r.file_size;
        });
    }

    if (files.empty()) {
        std::cerr << "no valid files to process" << std::endl << std::flush;
        return 1;
    }

    // processing
    double wtime_start = omp_get_wtime();
        std::vector<seen_t> results;
        std::vector<ReadsKmersStats> all_kmers_stats;

        RunInfo run_info;
        for (const auto& file_info : files) {
            std::cerr << " processing " << file_info.file_name << " size " << file_info.file_size << std::endl << std::flush;

            ReadsKmersStats kmers_stats;
            TimeStamp ts;
             
            // main procesing function. see processing.h/.cpp
            if(!process_file(file_info, run_info, kmers_stats, ts, results)) {
                std::cerr << "error processing file" << std::endl << std::flush;
            } else {
                std::cerr << file_info.file_name << " unique kmers: " << results.rbegin()->size() << std::endl << std::flush;
                std::cerr << kmers_stats;
                std::cerr << ts;
            }

            all_kmers_stats.push_back(kmers_stats);
        }

        if (results.empty()) {
            std::cerr << " no common kmers" << std::endl << std::flush;
            return 0;
        }

    // median sorting using heap
    double wtime_sort_start = omp_get_wtime();
        std::cerr << "median sorting" << std::endl << std::flush;
        auto cmp = [](const tuple_t& l, const tuple_t &r) {
            return std::get<1>(l) < std::get<1>(r);
        };
        std::priority_queue<tuple_t, std::vector<tuple_t>, decltype(cmp)> stat(cmp);
        auto last = *results.rbegin();
        for (seen_t::const_iterator it = last.begin(); it != last.end(); it++){
            std::vector<counter_t> tmp;
            for (auto r = results.begin(); r != results.end(); ++r) {
                const auto found = r->find(it->first);
                if (found != r->end()) {
                    tmp.push_back(found->second);
                }
            }
            std::vector<counter_t> md(tmp);
            std::nth_element(md.begin(), md.begin() + md.size()/2, md.end());

            stat.push(std::move(std::make_tuple(it->first, md[md.size()/2], std::move(tmp))));
        }

    // dumping final k-mers stats to std::cout
    double wtime_dump_start = omp_get_wtime();
        std::cerr << "dumping kmers" << std::endl << std::flush;
        std::vector<std::string> colnames;
        colnames.push_back("#kmer");
        for (auto &f : files)
            colnames.push_back(f.file_name);
        std::copy(colnames.begin(), colnames.end(), std::ostream_iterator<std::string>(std::cout, "\t"));
        std::cout << "median" << std::endl;

        while (!stat.empty()) {
           const auto &top = stat.top();

           char buf[KMER_LEN+1] = {0};
           from_u64mer(std::get<0>(top), buf, KMER_LEN); 

           const auto &med = std::get<1>(top);
           const auto &counts = std::get<2>(top);

           std::cout << buf << "\t";
           std::copy(counts.begin(), counts.end(), std::ostream_iterator<counter_t>(std::cout, "\t"));
           std::cout << med << std::endl;

           stat.pop();
    } 
    double wtime_end = omp_get_wtime();
    
    // times for benchmarking
    std::cerr << "#ALLT all " << wtime_end - wtime_start << std::endl << std::flush;
    std::cerr << "#ALLT process " << wtime_sort_start - wtime_start << std::endl << std::flush;
    std::cerr << "#ALLT sort " << wtime_dump_start - wtime_sort_start << std::endl << std::flush;
    std::cerr << "#ALLT dump " << wtime_end - wtime_dump_start << std::endl << std::flush;

    return 0;
}


// ext clone: git clone https://github.com/sparsehash/sparsehash-c11.git ext/sparsehash-c1
// build: g++ -std=c++14 -m64 -fopenmp  -O4 -I ext/sparsehash-c11 k31mare.cpp  processing.cpp io.cpp -o k31mare
// run: OMP_NUM_THREADS=4 /usr/bin/time -v k31mare data/*fastq 2>log > kmers

