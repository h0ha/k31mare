// io.h -- mmap magic and file info and runtime info classes

#ifndef _IO_H_
#define _IO_H_


#include <iostream>

#include <fcntl.h>
#include <omp.h>
#include <sys/mman.h> // mmap
#include <sys/stat.h> // fstat for file size
#include <unistd.h> // sysconf for offset

#include "config.h"


// various techincal information about run.
// number of threads. page sizes. default kmers stat storage size.
struct RunInfo {
    size_t page_size = 0;
    off_t offset_mask = 0;
    size_t threads = 0;
    size_t dict_size = DEFAULT_DICT_SIZE;

    RunInfo()
        : page_size(sysconf(_SC_PAGE_SIZE))
        {
            offset_mask = ~(page_size - 1);
            #pragma omp parallel default(none) shared(threads)
            #pragma omp master
            threads = omp_get_num_threads();
        }
};

// objects used for storing files names, sizes, fds
struct FileInfo {
    const char* file_name = NULL;
    size_t file_size = 0;
    int fd = -1;
    bool ok = false;

    // no default copy constructor, only mooving
    FileInfo(const FileInfo& fi) = delete;

    // constructor used for moving
    FileInfo(FileInfo&& fi)
        : file_name(std::move(fi.file_name))
        , file_size(std::move(fi.file_size))
        , fd(std::move(fi.fd))
        , ok(std::move(fi.ok))
        {
            fi.ok = false;
            fi.fd = -1;
            fi.file_size = 0;
            fi.file_name = NULL;
        }

    // assign operator used for moving
    FileInfo& operator= (FileInfo&& fi) {
        file_name = std::move(fi.file_name);
        file_size = std::move(fi.file_size);
        fd = std::move(fi.fd);
        ok = std::move(fi.ok);
       
        fi.ok = false;
        fi.fd = -1;
        fi.file_size = 0;
        fi.file_name = NULL;
    }
    
    FileInfo(const char* _file_name)
        : file_name(_file_name)
        {
            ok = false;
            if ((fd = open(file_name, O_RDONLY)) == -1) {
               std::cerr << "failed to open " << file_name
                         << std::endl << std::flush;
            } else {
                struct stat file_stats_buf;
                if (fstat(fd, &file_stats_buf) == -1) {
                    std::cerr << "failed to get filesize " << file_name
                              << std::endl << std::flush;
                    close(fd);
                    fd = -1;
                } else {
                    file_size = file_stats_buf.st_size; 
                    ok = true;
                } 
            }
        }

    ~FileInfo() {
        if (fd >= 0)
            close(fd);
        ok = false;
    }

    explicit operator bool() const { return ok; }
};


// objects used for storing mmap-io parameters, sizes, offsets, etc
struct BufInfo {
    bool mlocked = false;
    bool mmaped = false;
    bool ok = true;

    char *buf = NULL;
    size_t size = 0;
    size_t inbuf_size = 0;

    off_t pa_offset = 0;
    off_t inbuf_offset = 0;

    const RunInfo &run_info;
    const FileInfo &file_info;

    BufInfo(const RunInfo &_run_info, const FileInfo &_file_info)
        : run_info(_run_info)
        , file_info(_file_info)
        {
        }

    ~BufInfo() {
        if (mlocked)
            munlock(buf, size);
        mlocked = false;

        if (mmaped)
            munmap(buf, size);
        mmaped = false;
    }

    // mmap-io magic. mmaping, mlocking,  sizes / offsets handling. should be run
    void prepare_buf(off_t offset);

    // mmap-io magic. unmaping, unlocking, offsets handling. shoul be called at the end of iteration
    inline void free_buf(off_t &offset) {
        offset += inbuf_size;
        if (mlocked)
            munlock(buf, size);
        mlocked = false;

        if (mmaped)
            munmap(buf, size);
        mmaped = false;
    } 

    explicit operator bool() const { return ok; }
};



#endif // _IO_H_
