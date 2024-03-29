# CC = g++
# CFLAGS = -O1 -fno-tree-vectorize -Wall -Wextra -Wpedantic # -g -pg

CC = pgc++
CFLAGS = -Wall -Wextra -acc -fast -Minfo=all

FILES = read_hyperplanes lsh_single_table lsh_fast_creation lsh_fast_lookup \
		lsh_fast_creation_old

.PHONY: default all clean

default: lsh_fast_creation

all: $(FILES)

%: %.cpp
	$(CC) $(CFLAGS) $^ -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

lsh_single_table: lsh_single_table.o read_vector_file_funcs.o
	$(CC) $(CFLAGS) $^ -o $@

lsh_fast_creation: lsh_fast_creation.o read_vector_file_funcs.o
	$(CC) $(CFLAGS) $^ -o $@

lsh_fast_creation_old: lsh_fast_creation_old.o read_vector_file_funcs.o
	$(CC) $(CFLAGS) $^ -o $@

lsh_fast_lookup: lsh_fast_lookup.o read_vector_file_funcs.o
	$(CC) $(CFLAGS) $^ -o $@


run_lsh_single_table_siftsmall: lsh_single_table
	./lsh_single_table \
	../test_data/siftsmall/siftsmall_base.fvecs \
	../test_data/siftsmall/siftsmall_query.fvecs \
	../test_data/siftsmall/siftsmall_groundtruth.ivecs \
	../creating_hyperplanes/better_hyperplanes.dat

run_lsh_fast_creation_siftsmall: lsh_fast_creation
	./lsh_fast_creation \
	../test_data/siftsmall/siftsmall_base.fvecs \
	../test_data/siftsmall/siftsmall_query.fvecs \
	../test_data/siftsmall/siftsmall_groundtruth.ivecs \
	../creating_hyperplanes/better_hyperplanes.dat

run_lsh_fast_creation_siftsmall_old: lsh_fast_creation_old
	./lsh_fast_creation_old \
	../test_data/siftsmall/siftsmall_base.fvecs \
	../test_data/siftsmall/siftsmall_query.fvecs \
	../test_data/siftsmall/siftsmall_groundtruth.ivecs \
	../creating_hyperplanes/better_hyperplanes.dat

run_lsh_fast_lookup_siftsmall: lsh_fast_lookup
	./lsh_fast_lookup \
	../test_data/siftsmall/siftsmall_base.fvecs \
	../test_data/siftsmall/siftsmall_query.fvecs \
	../test_data/siftsmall/siftsmall_groundtruth.ivecs \
	../creating_hyperplanes/better_hyperplanes.dat

run_lsh_single_table_sift: lsh_single_table
	./lsh_single_table \
	../test_data/sift/sift_base.fvecs \
	../test_data/sift/sift_query.fvecs \
	../test_data/sift/sift_groundtruth.ivecs \
	../creating_hyperplanes/better_hyperplanes.dat

run_lsh_fast_creation_sift: lsh_fast_creation
	./lsh_fast_creation \
	../test_data/sift/sift_base.fvecs \
	../test_data/sift/sift_query.fvecs \
	../test_data/sift/sift_groundtruth.ivecs \
	../creating_hyperplanes/better_hyperplanes.dat

run_lsh_fast_creation_sift_old: lsh_fast_creation_old
	./lsh_fast_creation_old \
	../test_data/sift/sift_base.fvecs \
	../test_data/sift/sift_query.fvecs \
	../test_data/sift/sift_groundtruth.ivecs \
	../creating_hyperplanes/better_hyperplanes.dat

run_lsh_fast_lookup_sift: lsh_fast_lookup
	./lsh_fast_lookup \
	../test_data/sift/sift_base.fvecs \
	../test_data/sift/sift_query.fvecs \
	../test_data/sift/sift_groundtruth.ivecs \
	../creating_hyperplanes/better_hyperplanes.dat




clean:
	rm -f $(FILES) *.o
