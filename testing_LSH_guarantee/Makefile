# CC = g++
# CFLAGS = -O1 -fno-tree-vectorize -Wall -Wextra -Wpedantic # -g -pg

CC = pgc++
CFLAGS = -Wall -Wextra -acc -fast -Minfo=all
# CFLAGS = -Wall -Wextra -acc -fast -Minfo=all -lnvToolsExt

FILES = lsh_guarantee lsh_parallel_3 lsh_parallel_4

.PHONY: default all clean

default: lsh_guarantee

all: $(FILES)

%: %.cpp
	$(CC) $(CFLAGS) $^ -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

lsh_guarantee: lsh_guarantee.o read_file_funcs.o
	$(CC) $(CFLAGS) $^ -o $@

lsh_parallel_3: lsh_parallel_3.o read_file_funcs.o
	$(CC) $(CFLAGS) $^ -o $@

lsh_parallel_4: lsh_parallel_4.o read_file_funcs.o
	$(CC) $(CFLAGS) $^ -o $@


run_lsh_guarantee_siftsmall: lsh_guarantee
	./lsh_guarantee \
	../test_data/siftsmall/siftsmall_base.fvecs \
	../test_data/siftsmall/siftsmall_query.fvecs \
	../test_data/siftsmall/siftsmall_groundtruth.ivecs \
	../creating_hyperplanes/better_hyperplanes.dat

run_lsh_guarantee_sift: lsh_guarantee
	./lsh_guarantee \
	../test_data/sift/sift_base.fvecs \
	../test_data/sift/sift_query.fvecs \
	../test_data/sift/sift_groundtruth.ivecs \
	../creating_hyperplanes/better_hyperplanes.dat

profile_lsh_guarantee_sift_nvprof: lsh_guarantee
	nvprof --unified-memory-profiling off ./lsh_guarantee \
	../test_data/sift/sift_base.fvecs \
	../test_data/sift/sift_query.fvecs \
	../test_data/sift/sift_groundtruth.ivecs \
	../creating_hyperplanes/random_1000_zerosum_unit_hyperplanes.dat


run_lsh_parallel_3_siftsmall: lsh_parallel_3
	./lsh_parallel_3 \
	../test_data/siftsmall/siftsmall_base.fvecs \
	../test_data/siftsmall/siftsmall_query.fvecs \
	../test_data/siftsmall/siftsmall_groundtruth.ivecs \
	../creating_hyperplanes/better_hyperplanes.dat

run_lsh_parallel_3_sift: lsh_parallel_3
	./lsh_parallel_3 \
	../test_data/sift/sift_base.fvecs \
	../test_data/sift/sift_query.fvecs \
	../test_data/sift/sift_groundtruth.ivecs \
	../creating_hyperplanes/better_hyperplanes.dat

profile_lsh_parallel_3_sift_nvprof: lsh_parallel_3
	nvprof --unified-memory-profiling off ./lsh_parallel_3 \
	../test_data/sift/sift_base.fvecs \
	../test_data/sift/sift_query.fvecs \
	../test_data/sift/sift_groundtruth.ivecs \
	../creating_hyperplanes/random_1000_zerosum_unit_hyperplanes.dat


run_lsh_parallel_4_siftsmall: lsh_parallel_4
	./lsh_parallel_4 \
	../test_data/siftsmall/siftsmall_base.fvecs \
	../test_data/siftsmall/siftsmall_query.fvecs \
	../test_data/siftsmall/siftsmall_groundtruth.ivecs \
	../creating_hyperplanes/better_hyperplanes.dat

run_lsh_parallel_4_sift: lsh_parallel_4
	./lsh_parallel_4 \
	../test_data/sift/sift_base.fvecs \
	../test_data/sift/sift_query.fvecs \
	../test_data/sift/sift_groundtruth.ivecs \
	../creating_hyperplanes/better_hyperplanes.dat

profile_lsh_parallel_4_sift_nvprof: lsh_parallel_4
	nvprof --unified-memory-profiling off ./lsh_parallel_4 \
	../test_data/sift/sift_base.fvecs \
	../test_data/sift/sift_query.fvecs \
	../test_data/sift/sift_groundtruth.ivecs \
	../creating_hyperplanes/random_1000_zerosum_unit_hyperplanes.dat





clean:
	rm -f $(FILES) *.o
