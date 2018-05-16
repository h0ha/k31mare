== INTRO ==

This tools is written to process 31-mers.
All k-mers containig non-ACGT symbols are assumed to be invalid and discarded.
31-mers then packed into u64int_t with the 64th bit beeing used for error checking/flagging.
I don't have access to the icc and xeon phi thus no '#pragma vector' is used and no cache-tricks were used.  Only OpenMP and mmap io (input only) were used.


== BUILD ==

I'm building in 
    * ubuntu x86_64
with
    * gcc version 7.2.0 (Ubuntu 7.2.0-8ubuntu3.2) 

I'm using https://github.com/sparsehash/sparsehash-c11.git .
Should be already cloned in ext if not see below.

run the following commands:

# untar 
tar -zxf k31mare.tgz
cd k31mare

# clone sparce hash if not in ext
git clone https://github.com/sparsehash/sparsehash-c11.git ext/sparsehash-c11

# compile
g++ -std=c++14 -m64 -fopenmp  -O4 -I ext/sparsehash-c11 k31mare.cpp  processing.cpp io.cpp -o k31mare


== PREPARE ENVIRONMENT ==

Be sure either disable mlocking in config.h, or have enough 
'max locked memory       (kbytes, -l) 262144' ulimit.

i.e.
/etc/security/limits.conf:
user1            hard    memlock	 262144	
user1            soft    memlock	 262144	




== DATA == 

ftp://public-ftp.hmpdacc.org/Illumina/PHASEII/supragingival_plaque/SRS074976.tar.bz2
ftp://public-ftp.hmpdacc.org/Illumina/PHASEII/buccal_mucosa/SRS056461.tar.bz2

I've untared and moved everything into the single dir called `data`.


== RUNNIG ==

Hmm. So I have a 4 core laptop with 8Gb RAM and 16Gb SSD swap. And the idea was somehow to process 6 fastq files of 87G in total in reasonable time :)

Elapsed (wall clock) time (h:mm:ss or m:ss): 34:07.07
Maximum resident set size (kbytes): 7228016
Swaps: 0

This was done by

OMP_NUM_THREADS=4 /usr/bin/time -v \
    ~/projects/k31mare/k31mare --no-size-sort \
        data/SRS056461.denovo_duplicates_marked.trimmed.singleton.fastq \ 
        data/SRS074976.denovo_duplicates_marked.trimmed.singleton.fastq \
        data/SRS056461.denovo_duplicates_marked .trimmed.1.fastq \
        data/SRS056461.denovo_duplicates_marked.trimmed.2.fastq \
        data/SRS074976.denovo_duplicates_marked.trimmed.1.fastq \
        data/SRS074976.denovo_duplicates_marked.trimmed.2.fastq \
        2> log  > kmers

gzip log kmers

--no-size-sort leads to processing files in the specified order. And the 1st file (data/SRS056461.denovo_duplicates_marked.trimmed.singleton.fastq) has fewer kmvalid kmers than other files.

If You have enough memory just run:
OMP_NUM_THREADS=4 /usr/bin/time -v k31mare data/*fastq 2>log > kmers 

So I have 14298347 kmers with the following top-10:

#kmer	data/SRS056461.denovo_duplicates_marked.trimmed.singleton.fastq	data/SRS074976.denovo_duplicates_marked.trimmed.singleton.fastq	data/SRS056461.denovo_duplicates_marked.trimmed.1.fastq	data/SRS056461.denovo_duplicates_marked.trimmed.2.fastq	data/SRS074976.denovo_duplicates_marked.trimmed.1.fastq	data/SRS074976.denovo_duplicates_marked.trimmed.2.fastq	median
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG	364	3784	455	2803	9336	23296	3784
GGGTTGCGCTCGTTGCGGGACTTAACCCAAC	468	519	2022	1670	2339	2096	2022
GGTTGCGCTCGTTGCGGGACTTAACCCAACA	468	517	2022	1676	2347	2106	2022
TTACCCGACAAGGAATTTCGCTACCTTAGGA	440	573	2012	1541	2519	2180	2012
GTTGCGCTCGTTGCGGGACTTAACCCAACAT	470	513	2011	1660	2332	2068	2011
ACTTACCCGACAAGGAATTTCGCTACCTTAG	443	576	2006	1526	2527	2171	2006
CTTACCCGACAAGGAATTTCGCTACCTTAGG	442	571	2001	1526	2505	2154	2001
AACTTACCCGACAAGGAATTTCGCTACCTTA	442	577	1989	1528	2513	2166	1989
GGAACTTACCCGACAAGGAATTTCGCTACCT	437	580	1986	1529	2497	2161	1986
GAACTTACCCGACAAGGAATTTCGCTACCTT	436	574	1975	1523	2488	2150	1975


== BENCHMARKING ==

Two procedures were used.
1) Just looking at speed up.
The seed.fastq were generated contaning final kmers as reads. The file was put into the data dir and runs with various OMP_NUM_THREADS were performed.


zcat kmers.gz |
    tail -n +2 |
    awk '{id = "SRR00_"$8"_"NR; print "@"id; print $1; print "+"id; print "###############################" }' > data/seed.fastq

mkdir -p stats.seed
for thr in 1 2 3 4; do
    sync; sync
    OMP_NUM_THREADS=${thr} ~/projects/k31mare/k31mare data/*fastq 2>stats.seed/${thr}.log > stats.seed/${thr}.kmers
    sync; sync
    gzip -r stats.seed 
done


2) Scalability.
Reads from file were taken proportionally to the number of threads:
5% 10% 15% 20%.
From the top / no random reordering was used. 

Getting number of reads:
grep -F -c '@SRR' data/SRS0*fastq > data_read_counts

#data_read_counts 
# data/SRS056461.denovo_duplicates_marked.trimmed.1.fastq:90337316
# data/SRS056461.denovo_duplicates_marked.trimmed.2.fastq:90337316
# data/SRS056461.denovo_duplicates_marked.trimmed.singleton.fastq:22813617
# data/SRS074976.denovo_duplicates_marked.trimmed.1.fastq:83712382
# data/SRS074976.denovo_duplicates_marked.trimmed.2.fastq:83712382
# data/SRS074976.denovo_duplicates_marked.trimmed.singleton.fastq:20882055

or 

zcat log.sorted.gz |
    grep ^data/ -A 1 | grep total -B 1 |
    perl -pe 's/\s*$/ /g' | perl -pe 's/--/\n/g; s/$/\n/' |
    perl -pe 's/^\s+//' |
    cut -f 1,8 -d ' ' | perl -pe 's/ /:/' > data_read_counts_



mkdir -p stats.scal
for thr in 1 2 3 4; do
    perc=$((5 * ${thr}))
    mkdir -p data.scal
    rm -r data.scal
    mkdir -p data.scal 
    echo $perc
    for file in data/SRS*fastq; do
        total_reads=$(cat data_read_counts | grep -F $file | cut -f 2 -d ':')
        reads=$((${perc} * ${total_reads} / 100))   
        head -n $((4 * reads)) $file > data.scal/$(echo $file | cut -f 2 -d / )
    done

    sync; sync
    OMP_NUM_THREADS=${thr} ~/projects/k31mare/k31mare data.scal/*fastq 2>stats.scal/${thr}.log > stats.scal/${thr}.kmers
    sync; sync

    gzip -r stats.scal 
done


== PLOTTING == 

zgrep -F '#' stats.*.*/*log.gz |
    perl -pe 's,/,\t,g;
              s/:#(ALL)?T\s+/\tT$1\t/;
              s/^stats.//;
              s/\./\t/ ;
              s/.log.gz//;
              s/\s+([\d\.\-e]+)\s*$/\t$1\n/;
              s/\s*:\s*\t/\t/' > times.tsv
    

# R
library(ggplot2)
library(data.table)

raw = read.delim("times.tsv",sep = "\t", as.is = T, header = F);
raw = data.table(raw)
colnames(raw) = c("bench", "rep", "threads", "type", "action", "time")

data = raw[,.(time = sum(time)), by = c("bench", "rep", "threads", "type", "action")]

data_m_pre_1 = merge(data[action != "hash merging"], data[threads == 1 & action != "hash merging"], all = T,  by = c("bench", "rep", "type", "action"))
data_m_pre_2 = merge(data[action == "hash merging"], data[threads == 2 & action == "hash merging"], all = T,  by = c("bench", "rep", "type", "action"))

data_m = rbind(data_m_pre_1, data_m_pre_2) 

data_scal_m = data_m[bench == "scal"]
data_scal_m$scal =  data_scal_m$time.x / (data_scal_m$time.y * data_scal_m$threads.x) 

colnames(data_scal_m)
[1] "bench"     "rep"       "type"      "action"    "threads.x" "time.x"   
[7] "threads.y" "time.y"    "scal"     

ggplot(data_scal_m[type=="TALL"], aes(y = scal, x = threads.x, color = action, shape = action)) + geom_point(size = 4, alpha = 0.7, position = position_jitter(width = 0.05, height = 0)) + geom_smooth(alpha = 0.2, method=lm) + theme_bw() + labs(title = "Scalability (with median sorting and output)", y = "Sclability", x = "Threads")
ggsave("scalability.everything.png")

ggplot(data_scal_m[type=="T"], aes(y = scal, x = threads.x, color = action, shape = action)) + geom_point(size = 4, alpha = 0.7, position = position_jitter(width = 0.05, height = 0)) + geom_smooth(alpha = 0.2, method=lm) + theme_bw() + scale_shape_manual(values =c(15:20, 21:25)) + labs(title = "Scalability (just processing, no median sorting and output)", y = "Scalability", x = "Threads")
ggsave("scalability.processing.png")

ggplot(data_scal_m[type=="T" & threads.x > 1], aes(y = scal, x = threads.x, color = action, shape = action)) + geom_point(size = 4, alpha = 0.7, position = position_jitter(width = 0.05, height = 0)) + geom_smooth(alpha = 0.2, method=lm) + theme_bw() + scale_shape_manual(values =c(15:20, 21:25)) + labs(title = "Scalability (just processing, no median sorting and output, threads > 1)", y = "Scalability", x = "Threads")
ggsave("scalability.processing.threads.gt1.png")



data_seed_m = data_m[bench == "seed"]
data_seed_m$accel =  data_seed_m$time.y / data_seed_m$time.x 

data_seed_m = data_seed_m[!(action == "hash merging" & threads.x == 1)] 

ggplot(data_seed_m[type=="TALL"], aes(y = accel, x = threads.x, color = action, shape = action)) + geom_point(size = 4, alpha = 0.7, position = position_jitter(width = 0.05, height = 0)) + geom_smooth(alpha = 0.2, method=lm) + theme_bw() + labs(title = "Acceleration (with median sorting and output)", y = "Acceleration", x = "Threads") + geom_abline(slope = 1, color = "red")
ggsave("acceleration.everything.png")

ggplot(data_seed_m[type=="T"], aes(y = accel, x = threads.x, color = action, shape = action)) + geom_point(size = 4, alpha = 0.7, position = position_jitter(width = 0.05, height = 0)) + geom_smooth(alpha = 0.05, method=lm) + theme_bw() + scale_shape_manual(values =c(15:20, 21:25)) + labs(title = "Acceleration (just processing, no median sorting and output)", y = "Acceleration", x = "Threads") + geom_abline(slope = 1, color = "red")
ggsave("acceleration.processing.png")



