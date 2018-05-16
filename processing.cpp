// processing.cpp -- main kmers processing functions

#include <string.h> // for memchr

#include "processing.h"
#include "u64mer.h"
#include <algorithm>


// main file processing function
// OpenMP threading and mmap-io is used in it
bool process_file(const FileInfo &file_info, const RunInfo &run_info,
         ReadsKmersStats &kmers_stats, TimeStamp &ts,
         std::vector<seen_t>& results) {

    ts.start_epoch();

    size_t dict_size = run_info.dict_size;

    //trying to find out if it's possible to use results from the previous run as filter
    const seen_t *filter = NULL;
    if (!results.empty()) {
        filter = &results[results.size()-1];
        dict_size = std::min(dict_size, filter->size());
        std::cerr << "using filter " << filter->size() << " large" << std::endl << std::flush;
    }


    // initialize kmers storages. separate for each thread
    std::vector<seen_t> seen(run_info.threads, seen_t(dict_size));
    #pragma omp parallel default(none) shared(seen, dict_size)
    {
        size_t thrid = omp_get_thread_num();
        seen[thrid].set_deleted_key(1UL<<63);
        seen[thrid].resize(dict_size);
        //seen[thrid].min_load_factor(0.0);
    }

    ts.start_iterations();

    size_t iter = 0;
    off_t offset = 0;
    BufInfo bif(run_info, file_info);

    // reads, kmers stats special reduction
    #pragma omp declare reduction(sumReadsKmers : ReadsKmersStats : omp_out += omp_in)


    // main processing loop
    #pragma omp parallel default(none) \
        shared(run_info, file_info, bif, ts)\
        shared(offset, iter, seen, filter, std::cerr) \
        reduction(sumReadsKmers: kmers_stats)
    while (file_info.file_size && offset < file_info.file_size) {
        size_t threads = omp_get_num_threads();
        size_t thrid = omp_get_thread_num();

        #pragma omp master
        {
            ts.start_iteration();
            ts.start_io();

            // map io magic, mmap, lock, sizes/offsets calculation
            bif.prepare_buf(offset);

            ts.end_io();
        }

        #pragma omp barrier
        if (!bif)
            break;

        #pragma omp master
        ts.start_processing();

        // find out 'true'  borders of reads chunk
        // form proper/valid reads chunck
        char *first = bif.buf + bif.inbuf_offset;
        char *last  = bif.buf + bif.inbuf_offset + bif.inbuf_size - 1;

        if (bif.inbuf_size <= MIN_DATA_SIZE_IN_PAGES * run_info.page_size) {
            if (thrid != 0) {
                 first = NULL;
                 last = NULL;
            }
        } else {
            size_t chunk_size = bif.inbuf_size / threads;
            if (thrid != 0) {
                first = 1 + find_full_read_end(bif.buf + bif.inbuf_offset, chunk_size * thrid);
            }
            if (thrid != (threads-1)) {
                last = find_full_read_end(bif.buf + bif.inbuf_offset, chunk_size * (thrid+1));
            }
        }

        // process chunck of reads
        if (first != last) {
            fill_stats(first, last, seen[thrid], kmers_stats, filter);
        }

        #pragma omp barrier

        #pragma omp master
        {
            ts.end_processing();
            ts.start_io();
            
            // map-io magic, unmap, unlock, increment offset 
            bif.free_buf(offset);
            
            ++iter;

            // dump master thread biased cryptic info, aka logging / keep-alive 'signals'
            if (iter % PRINT_INFO_EACH == 0) {
                std::cerr << kmers_stats;
                std::cerr << "master unique kmers: " << seen[thrid].size() << std::endl << std::flush;
            }

            ts.end_io();
            ts.end_iteration();
        }

        #pragma omp barrier
    }

    ts.end_iterations();

    if (!bif)
        return false;

    #pragma omp parallel
    {
        size_t thrid = omp_get_thread_num();
        #pragma omp critical
        std::cerr << thrid << " unique kmers: " << seen[thrid].size() << std::endl << std::flush;
    }


    // merging k-mers stats from thread's storages into the master/0 one
    ts.start_merge();

    std::cerr << "merging kmer maps" << std::endl << std::flush;

    for (auto thrid = 1; thrid < run_info.threads; thrid++) {
        std::cerr << "merging " << thrid << std::endl << std::flush;
        seen[thrid].resize(seen[thrid].size());

        for (seen_t::const_iterator it = seen[thrid].begin(); it != seen[thrid].end(); ++it)
            seen[0][it->first] += it->second;
        seen[thrid].clear();
        seen[thrid].resize(0);

        seen_t tmpl;
        tmpl.swap(seen[thrid]);
    }
    seen[0].resize(seen[0].size());

    ts.end_merge();


    // trying to remove unused kmers from the previous resuls
    ts.start_filter();

    std::cerr << "filtering prev kmer map" << std::endl << std::flush;
    if (!seen[0].empty() && !results.empty() ) {
        const auto& cur = seen[0];

        auto i = 0;
        for (auto prev = results.begin(); prev != results.end(); ++prev, ++i) {
            std::cerr << "filtering for " << i << std::endl << std::flush;
            prev->set_deleted_key(1UL<<63);
            for (seen_t::iterator it = prev->begin(); it != prev->end(); ++it) {
                if (cur.find(it->first) == cur.end()) {
                    prev->erase(it->first);
                }
            }
            prev->resize(prev->size());
        }
    }

    ts.end_filter();

    // use merged k-mers as the result
    if (!seen[0].empty()) {
        results.push_back(seen_t());
        results.rbegin()->swap(seen[0]);
    }
    seen.clear();

    ts.end_epoch();

    return true;
}



// fill k-mers into stat from the well formed chunk of reads
// assumed that start points to the proper read start and reads ids are starting from the same prefix of READID_PFX_LEN, see config.h
void fill_stats(char *start, char *end, seen_t &seen, ReadsKmersStats &kmers_stats, const seen_t *filter, char eol) {
    char *start_initial = start;

    auto line = 0;
    char *eol_ptr = (char*)memchr(start, eol, end - start + 1);
    while (start < end && eol_ptr != NULL) {

        #ifdef CHECK_READ_START
        if (memcmp(start, start_initial, READID_PFX_LEN) == 0)
            line = 0;
        #endif // CHECK_READ_START

        size_t len = eol_ptr + 1 - start;

        if (line == 0) ++kmers_stats.total_reads;

        // pack into u64mer and store
        if (line == 1 && len >= KMER_LEN) {
            for (off_t i = 0; i < len - KMER_LEN; i++) {
                ++kmers_stats.total_kmers;

                uint64_t res = u64mer(start + i);
                if (test_u64mer(res)) {
                    if (filter == NULL || filter->find(res) != filter->end()) {
                        ++seen[res];
                    }
                } else {
                    ++kmers_stats.bad_kmers;
                }
            }
        }
        
        line = (line + 1) & (4-1);
        start = eol_ptr+1;
        eol_ptr = (char*)memchr(start, eol, end - start);
    }
}


// util used to find the farest good read border
// assumed thar full read consists of 4 lines ending with EOL (see config.h)
/* what if the chunk / buffers ends with the improper read's part. ie.:
    [\n] or start
    @srr\n
    acgt\n
    +srr\n
    qual\n
    @sr
*/
char *find_full_read_end(char *buf, size_t size, char eolc){
    char *eop_last = NULL;

    // searching from the end of the buffer. looking only for the last 5 line or less.
    for (auto lines_back = 0; lines_back < 5; lines_back++) {
        char *eop = (char*) memrchr(buf, eolc, size);
        if (lines_back == 0) // store the last(first from the buffer end) eol position
            eop_last = eop;
        size = eop - buf;

        // if we have found start of the read 4 lines before
        // and read ends with the last eol position, the read was well-formed then
        // return its end
        if (lines_back == 4) {
            if(eop == NULL || memcmp(buf, eop+1, READID_PFX_LEN) == 0)
                return eop_last;
        }

        if (eop == NULL)
            return NULL;

        // return if we've found the start of the read
        if (memcmp(buf, eop+1, READID_PFX_LEN) == 0)
            return eop;
    }
    
    return NULL;
}
