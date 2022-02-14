CC = g++
# CFLAGS = -Wall -Wextra -Wpedantic
# CFLAGS = -O3 -Wall -Wextra -Wpedantic
CFLAGS = -O1 -fno-tree-vectorize -Wall -Wextra -Wpedantic -g -pg
FILES = lsh_serial lsh_serial_random_points mlsh_serial mlsh_serial2 \
	read_hyperplane_tables lsh_serial_loadHP mlsh_serial_loadHP

.PHONY: default all clean profile_lsh

default: all

all: $(FILES)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

%: %.cpp
	$(CC) $(CFLAGS) $^ -o $@

lsh_serial: lsh_serial.o read_file_funcs.o
	$(CC) $(CFLAGS) $^ -o $@

lsh_serial_loadHP: lsh_serial_loadHP.o read_file_funcs.o
	$(CC) $(CFLAGS) $^ -o $@

mlsh_serial_loadHP: mlsh_serial_loadHP.o read_file_funcs.o
	$(CC) $(CFLAGS) $^ -o $@


mlsh_serial: mlsh_serial.o read_file_funcs.o
	$(CC) $(CFLAGS) $^ -o $@

mlsh_serial2: mlsh_serial2.o read_file_funcs.o
	$(CC) $(CFLAGS) $^ -o $@

# profile_lsh:
# 	gprof lsh_serial gmon.out > gprof_results.txt

run_lsh_serial_siftsmall: lsh_serial
	./lsh_serial \
	../test_data/siftsmall/siftsmall_base.fvecs \
	../test_data/siftsmall/siftsmall_query.fvecs \
	../test_data/siftsmall/siftsmall_groundtruth.ivecs

run_lsh_serial_sift: lsh_serial
	./lsh_serial \
	../test_data/sift/sift_base.fvecs \
	../test_data/sift/sift_query.fvecs \
	../test_data/sift/sift_groundtruth.ivecs

run_mlsh_serial_siftsmall: mlsh_serial
	./mlsh_serial \
	../test_data/siftsmall/siftsmall_base.fvecs \
	../test_data/siftsmall/siftsmall_query.fvecs \
	../test_data/siftsmall/siftsmall_groundtruth.ivecs

run_mlsh_serial_sift: mlsh_serial
	./mlsh_serial \
	../test_data/sift/sift_base.fvecs \
	../test_data/sift/sift_query.fvecs \
	../test_data/sift/sift_groundtruth.ivecs

run_mlsh_serial_siftsmall2: mlsh_serial2
	./mlsh_serial2 \
	../test_data/siftsmall/siftsmall_base.fvecs \
	../test_data/siftsmall/siftsmall_query.fvecs \
	../test_data/siftsmall/siftsmall_groundtruth.ivecs

run_mlsh_serial_sift2: mlsh_serial2
	./mlsh_serial2 \
	../test_data/sift/sift_base.fvecs \
	../test_data/sift/sift_query.fvecs \
	../test_data/sift/sift_groundtruth.ivecs


run_lsh_serial_loadHP_siftsmall: lsh_serial_loadHP
	./lsh_serial_loadHP \
	../test_data/siftsmall/siftsmall_base.fvecs \
	../test_data/siftsmall/siftsmall_query.fvecs \
	../test_data/siftsmall/siftsmall_groundtruth.ivecs

run_lsh_serial_loadHP_sift: lsh_serial_loadHP
	./lsh_serial_loadHP \
	../test_data/sift/sift_base.fvecs \
	../test_data/sift/sift_query.fvecs \
	../test_data/sift/sift_groundtruth.ivecs

run_mlsh_serial_loadHP_siftsmall: mlsh_serial_loadHP
	./mlsh_serial_loadHP \
	../test_data/siftsmall/siftsmall_base.fvecs \
	../test_data/siftsmall/siftsmall_query.fvecs \
	../test_data/siftsmall/siftsmall_groundtruth.ivecs

run_mlsh_serial_loadHP_sift: mlsh_serial_loadHP
	./mlsh_serial_loadHP \
	../test_data/sift/sift_base.fvecs \
	../test_data/sift/sift_query.fvecs \
	../test_data/sift/sift_groundtruth.ivecs


clean:
	rm -f $(FILES) *.o
