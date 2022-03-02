# CC = g++
# # CFLAGS = -O3 -Wall -Wextra -Wpedantic
# CFLAGS = -O1 -fno-tree-vectorize -Wall -Wextra -Wpedantic -g -pg

CC = pgc++
CFLAGS = -Wall -Wextra -Wpedantic -acc -fast -Minfo=all
# CFLAGS = -Wall -Wextra -Wpedantic -Minfo=all -g -pg

FILES = lsh_serial lsh_parallel_3 lsh_parallel_4 lsh_parallel_5 lsh_parallel_6 lsh_parallel_7

.PHONY: default all clean profile_lsh

default: lsh_parallel_7

all: $(FILES)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

%: %.cpp
	$(CC) $(CFLAGS) $^ -o $@

lsh_serial: lsh_serial.o read_file_funcs.o
	$(CC) $(CFLAGS) $^ -o $@

lsh_parallel_3: lsh_parallel_3.o read_file_funcs.o
	$(CC) $(CFLAGS) $^ -o $@

lsh_parallel_4: lsh_parallel_4.o read_file_funcs.o
	$(CC) $(CFLAGS) $^ -o $@

lsh_parallel_5: lsh_parallel_5.o read_file_funcs.o
	$(CC) $(CFLAGS) $^ -o $@

lsh_parallel_6: lsh_parallel_6.o read_file_funcs.o
	$(CC) $(CFLAGS) $^ -o $@

lsh_parallel_7: lsh_parallel_7.o read_file_funcs.o
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

run_lsh_parallel_3_siftsmall: lsh_parallel_3
	./lsh_parallel_3 \
	../test_data/siftsmall/siftsmall_base.fvecs \
	../test_data/siftsmall/siftsmall_query.fvecs \
	../test_data/siftsmall/siftsmall_groundtruth.ivecs \
	../creating_hyperplanes/random_1000_zerosum_unit_hyperplanes.dat

run_lsh_parallel_3_sift: lsh_parallel_3
	./lsh_parallel_3 \
	../test_data/sift/sift_base.fvecs \
	../test_data/sift/sift_query.fvecs \
	../test_data/sift/sift_groundtruth.ivecs \
	../creating_hyperplanes/random_1000_zerosum_unit_hyperplanes.dat

profile_lsh_parallel_3_siftsmall: lsh_parallel_3
	nvprof --unified-memory-profiling off ./lsh_parallel_3 \
	../test_data/siftsmall/siftsmall_base.fvecs \
	../test_data/siftsmall/siftsmall_query.fvecs \
	../test_data/siftsmall/siftsmall_groundtruth.ivecs \
	../creating_hyperplanes/random_1000_zerosum_unit_hyperplanes.dat

profile_lsh_parallel_3_sift: lsh_parallel_3
	nvprof --unified-memory-profiling off ./lsh_parallel_3 \
	../test_data/sift/sift_base.fvecs \
	../test_data/sift/sift_query.fvecs \
	../test_data/sift/sift_groundtruth.ivecs \
	../creating_hyperplanes/random_1000_zerosum_unit_hyperplanes.dat


run_lsh_parallel_4_sift: lsh_parallel_4
	./lsh_parallel_4 \
	../test_data/sift/sift_base.fvecs \
	../test_data/sift/sift_query.fvecs \
	../test_data/sift/sift_groundtruth.ivecs \
	../creating_hyperplanes/random_1000_zerosum_unit_hyperplanes.dat

profile_lsh_parallel_4_sift: lsh_parallel_4
	nvprof --unified-memory-profiling off ./lsh_parallel_4 \
	../test_data/sift/sift_base.fvecs \
	../test_data/sift/sift_query.fvecs \
	../test_data/sift/sift_groundtruth.ivecs \
	../creating_hyperplanes/random_1000_zerosum_unit_hyperplanes.dat


run_lsh_parallel_5_sift: lsh_parallel_5
	./lsh_parallel_5 \
	../test_data/sift/sift_base.fvecs \
	../test_data/sift/sift_query.fvecs \
	../test_data/sift/sift_groundtruth.ivecs \
	../creating_hyperplanes/random_1000_zerosum_unit_hyperplanes.dat

profile_lsh_parallel_5_sift: lsh_parallel_5
	nvprof --unified-memory-profiling off ./lsh_parallel_5 \
	../test_data/sift/sift_base.fvecs \
	../test_data/sift/sift_query.fvecs \
	../test_data/sift/sift_groundtruth.ivecs \
	../creating_hyperplanes/random_1000_zerosum_unit_hyperplanes.dat


run_lsh_parallel_6_sift: lsh_parallel_6
	./lsh_parallel_6 \
	../test_data/sift/sift_base.fvecs \
	../test_data/sift/sift_query.fvecs \
	../test_data/sift/sift_groundtruth.ivecs \
	../creating_hyperplanes/random_1000_zerosum_unit_hyperplanes.dat

profile_lsh_parallel_6_sift: lsh_parallel_6
	nvprof --unified-memory-profiling off ./lsh_parallel_6 \
	../test_data/sift/sift_base.fvecs \
	../test_data/sift/sift_query.fvecs \
	../test_data/sift/sift_groundtruth.ivecs \
	../creating_hyperplanes/random_1000_zerosum_unit_hyperplanes.dat

# profile_lsh_parallel_6_sift_ncu: lsh_parallel_6
# 	sudo /opt/nvidia/hpc_sdk/Linux_x86_64/21.9/compilers/bin/ncu -o profile \
# 	./lsh_parallel_6 \
# 	../test_data/sift/sift_base.fvecs \
# 	../test_data/sift/sift_query.fvecs \
# 	../test_data/sift/sift_groundtruth.ivecs \
# 	../creating_hyperplanes/random_1000_zerosum_unit_hyperplanes.dat

# profile_lsh_parallel_6_sift_nsys: lsh_parallel_6
# 	sudo /opt/nvidia/hpc_sdk/Linux_x86_64/21.9/compilers/bin/nsys \
# 	./lsh_parallel_6 \
# 	../test_data/sift/sift_base.fvecs \
# 	../test_data/sift/sift_query.fvecs \
# 	../test_data/sift/sift_groundtruth.ivecs \
# 	../creating_hyperplanes/random_1000_zerosum_unit_hyperplanes.dat

run_lsh_parallel_7_sift: lsh_parallel_7
	./lsh_parallel_7 \
	../test_data/sift/sift_base.fvecs \
	../test_data/sift/sift_query.fvecs \
	../test_data/sift/sift_groundtruth.ivecs \
	../creating_hyperplanes/random_1000_zerosum_unit_hyperplanes.dat

profile_lsh_parallel_7_sift: lsh_parallel_7
	nvprof --unified-memory-profiling off ./lsh_parallel_7 \
	../test_data/sift/sift_base.fvecs \
	../test_data/sift/sift_query.fvecs \
	../test_data/sift/sift_groundtruth.ivecs \
	../creating_hyperplanes/random_1000_zerosum_unit_hyperplanes.dat


# ../test_data/sift/sift_base.fvecs ../test_data/sift/sift_query.fvecs ../test_data/sift/sift_groundtruth.ivecs ../creating_hyperplanes/random_1000_zerosum_unit_hyperplanes.dat
# ../test_data/siftsmall/siftsmall_base.fvecs ../test_data/siftsmall/siftsmall_query.fvecs ../test_data/siftsmall/siftsmall_groundtruth.ivecs ../creating_hyperplanes/random_1000_zerosum_unit_hyperplanes.dat


clean:
	rm -f $(FILES) *.o
