// io.cpp -- mmap magic. BufInfo::prepare_buf code
#include "io.h"
#include "processing.h"


void BufInfo::prepare_buf(off_t offset) {
    ok = true;

    // offset for mmap shoul be paged aligned (pa_)
    pa_offset = offset & run_info.offset_mask;
    inbuf_offset = offset & ~run_info.offset_mask;

    size = BUFFER_SIZE_IN_PAGES * run_info.page_size; 

    // size can't be larger the rest of the file
    if (pa_offset + size > file_info.file_size)
        size = file_info.file_size - pa_offset;

    // mmap
    buf = (char*)mmap(NULL, size,
        PROT_READ, MAP_PRIVATE | MAP_POPULATE,
        file_info.fd, pa_offset);

    if (buf == MAP_FAILED) {
        std::cerr << "mmap failed for " << file_info.file_name
                  << " offset " << offset << " size " << size << std::endl << std::flush; 

        ok = false;
        return;
    }
    mmaped = true;

    #ifdef USE_MLOCK
        // mlocking is safer for multithreading, 
        // but be sure to have enough ulimits for this
        if (mlock(buf, size) != 0) {
            std::cerr << "mlock failed for " << file_info.file_name
                  << " offset " << offset << " size " << size << std::endl << std::flush; 
            std::cerr << "try to set ulimit larger then " << size << " bytes" << std::endl << std::flush;

            if (mmaped)
                munmap(buf, size);
            mmaped = false;

            ok = false;
            return;
        }
        mlocked = true;
    #else
        mlocked = false;
    #endif
    

    // find the farest proper read's ending
    char *eol_ptr = find_full_read_end(buf + inbuf_offset, size-inbuf_offset);
    // will die even if the only one last is corrupted
    if (eol_ptr == NULL)  {
        std::cerr << "wrong format: need eol to process " << std::endl << std::flush;

        if (mlocked)
            munlock(buf, size);
        if (mmaped)
            munmap(buf, size);

        ok = false;
        return;
    }

    inbuf_size = eol_ptr+1 - (buf + inbuf_offset);
}

