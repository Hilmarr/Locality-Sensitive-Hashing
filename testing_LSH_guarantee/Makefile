# CC = g++
# CFLAGS = -O1 -fno-tree-vectorize -Wall -Wextra -Wpedantic # -g -pg

CC = pgc++
CFLAGS = -Wall -Wextra -acc -fast -Minfo=all

FILES = lsh_guarantee

.PHONY: default all clean

default: lsh_guarantee

all: $(FILES)

%: %.cpp
	$(CC) $(CFLAGS) $^ -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

lsh_guarantee: lsh_guarantee.o read_file_funcs.o
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



clean:
	rm -f $(FILES) *.o
