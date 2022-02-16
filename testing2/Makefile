# CC = g++
# # CFLAGS = -O3 -Wall -Wextra -Wpedantic
# CFLAGS = -O1 -fno-tree-vectorize -Wall -Wextra -Wpedantic -g -pg

CC = pgc++
CFLAGS = -Wall -Wextra -Wpedantic -acc -fast -Minfo=all
# CFLAGS = -Wall -Wextra -Wpedantic -Minfo=all -g -pg

FILES = lsh_serial lsh_parallel_1

.PHONY: default all clean profile_lsh

default: lsh_parallel_1

all: $(FILES)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

%: %.cpp
	$(CC) $(CFLAGS) $^ -o $@

lsh_serial: lsh_serial.o read_file_funcs.o
	$(CC) $(CFLAGS) $^ -o $@

lsh_parallel_1: lsh_parallel_1.o read_file_funcs.o
	$(CC) $(CFLAGS) $^ -o $@

profile_lsh_serial:
	gprof lsh_serial gmon.out > gprof_results.txt

run_lsh_serial_siftsmall: lsh_serial
	./lsh_serial \
	../test_data/siftsmall/siftsmall_base.fvecs \
	../test_data/siftsmall/siftsmall_query.fvecs \
	../test_data/siftsmall/siftsmall_groundtruth.ivecs \
	../creating_hyperplanes/random_1000_zerosum_unit_hyperplanes.dat

run_lsh_serial_sift: lsh_serial
	./lsh_serial \
	../test_data/sift/sift_base.fvecs \
	../test_data/sift/sift_query.fvecs \
	../test_data/sift/sift_groundtruth.ivecs \
	../creating_hyperplanes/random_1000_zerosum_unit_hyperplanes.dat

run_lsh_parallel_1_siftsmall: lsh_parallel_1
	./lsh_parallel_1 \
	../test_data/siftsmall/siftsmall_base.fvecs \
	../test_data/siftsmall/siftsmall_query.fvecs \
	../test_data/siftsmall/siftsmall_groundtruth.ivecs \
	../creating_hyperplanes/random_1000_zerosum_unit_hyperplanes.dat

run_lsh_parallel_1_sift: lsh_parallel_1
	./lsh_parallel_1 \
	../test_data/sift/sift_base.fvecs \
	../test_data/sift/sift_query.fvecs \
	../test_data/sift/sift_groundtruth.ivecs \
	../creating_hyperplanes/random_1000_zerosum_unit_hyperplanes.dat

clean:
	rm -f $(FILES) *.o
