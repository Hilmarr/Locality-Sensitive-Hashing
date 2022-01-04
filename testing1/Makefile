CC = g++
# CFLAGS = -Wall -Wextra -Wpedantic
# CFLAGS = -O3 -Wall -Wextra -Wpedantic
CFLAGS = -O1 -fno-tree-vectorize -Wall -Wextra -Wpedantic # -g -pg
FILES = lsh_serial

.PHONY: default all clean profile_lsh

default: lsh_serial

all: $(FILES)

%: %.cpp
	$(CC) $(CFLAGS) $^ -o $@

# profile_lsh:
# 	gprof lsh_serial gmon.out > gprof_results.txt

run_lsh_serial_siftsmall: lsh_serial
	./lsh_serial \
	../test_data/siftsmall/siftsmall_base.fvecs \
	../test_data/siftsmall/siftsmall_query.fvecs \
	../test_data/siftsmall/siftsmall_groundtruth.ivecs

clean:
	rm -f $(FILES) *.o
