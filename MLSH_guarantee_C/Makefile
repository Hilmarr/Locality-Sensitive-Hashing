CC = g++
CFLAGS = -O1 -fno-tree-vectorize -Wall -Wextra -Wpedantic # -g -pg
FILES = read_hyperplanes lsh_single_table

.PHONY: default all clean

default: lsh_single_table

all: $(FILES)

%: %.cpp
	$(CC) $(CFLAGS) $^ -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

lsh_single_table: lsh_single_table.o read_vector_file_funcs.o
	$(CC) $(CFLAGS) $^ -o $@

run_lsh_single_table_siftsmall: lsh_single_table
	./lsh_single_table \
	../test_data/siftsmall/siftsmall_base.fvecs \
	../test_data/siftsmall/siftsmall_query.fvecs \
	../test_data/siftsmall/siftsmall_groundtruth.ivecs \
	../hyperplanes.dat



clean:
	rm -f $(FILES) *.o
