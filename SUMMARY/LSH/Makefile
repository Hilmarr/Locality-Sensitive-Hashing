
CC = pgc++
CFLAGS = -Wall -Wextra -acc -fast -Minfo=all -lnvToolsExt

FILES = lsh_multi_table lsh_single_table

.PHONY: default all clean profile_lsh

default: all

all: $(FILES)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

%: %.cpp
	$(CC) $(CFLAGS) $^ -o $@


lsh_multi_table: lsh_multi_table.o read_file_funcs.o
	$(CC) $(CFLAGS) $^ -o $@

lsh_single_table: lsh_single_table.o read_file_funcs.o
	$(CC) $(CFLAGS) $^ -o $@



run_lsh_multi_table_sift: lsh_multi_table
	./lsh_multi_table \
	../../test_data/sift/sift_base.fvecs \
	../../test_data/sift/sift_query.fvecs \
	../../test_data/sift/sift_groundtruth.ivecs \
	../Hyperplanes/random_1000_zerosum_unit_hplanes

profile_lsh_multi_table_sift: lsh_multi_table
	nvprof --unified-memory-profiling off ./lsh_multi_table \
	../../test_data/sift/sift_base.fvecs \
	../../test_data/sift/sift_query.fvecs \
	../../test_data/sift/sift_groundtruth.ivecs \
	../Hyperplanes/random_1000_zerosum_unit_hplanes

run_lsh_multi_table_siftsmall: lsh_multi_table
	./lsh_multi_table \
	../../test_data/siftsmall/siftsmall_base.fvecs \
	../../test_data/siftsmall/siftsmall_query.fvecs \
	../../test_data/siftsmall/siftsmall_groundtruth.ivecs \
	../Hyperplanes/random_1000_zerosum_unit_hplanes

profile_lsh_multi_table_siftsmall: lsh_multi_table
	nvprof --unified-memory-profiling off ./lsh_multi_table \
	../../test_data/siftsmall/siftsmall_base.fvecs \
	../../test_data/siftsmall/siftsmall_query.fvecs \
	../../test_data/siftsmall/siftsmall_groundtruth.ivecs \
	../Hyperplanes/random_1000_zerosum_unit_hplanes


run_lsh_single_table_sift: lsh_single_table
	./lsh_single_table \
	../../test_data/sift/sift_base.fvecs \
	../../test_data/sift/sift_query.fvecs \
	../../test_data/sift/sift_groundtruth.ivecs \
	../Hyperplanes/better_hplanes

profile_lsh_single_table_sift: lsh_single_table
	nvprof --unified-memory-profiling off ./lsh_single_table \
	../../test_data/sift/sift_base.fvecs \
	../../test_data/sift/sift_query.fvecs \
	../../test_data/sift/sift_groundtruth.ivecs \
	../Hyperplanes/better_hplanes

run_lsh_single_table_siftsmall: lsh_single_table
	./lsh_single_table \
	../../test_data/siftsmall/siftsmall_base.fvecs \
	../../test_data/siftsmall/siftsmall_query.fvecs \
	../../test_data/siftsmall/siftsmall_groundtruth.ivecs \
	../Hyperplanes/better_hplanes

profile_lsh_single_table_siftsmall: lsh_single_table
	nvprof --unified-memory-profiling off ./lsh_single_table \
	../../test_data/siftsmall/siftsmall_base.fvecs \
	../../test_data/siftsmall/siftsmall_query.fvecs \
	../../test_data/siftsmall/siftsmall_groundtruth.ivecs \
	../Hyperplanes/better_hplanes



clean:
	rm -f $(FILES) *.o
