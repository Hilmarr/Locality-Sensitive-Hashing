CC = g++
CFLAGS = -O1 -fno-tree-vectorize -Wall -Wextra -Wpedantic
FILES = read_data func_read_data

.PHONY: default all clean run_read_data

default: func_read_data

all: $(FILES)

%: %.cpp
	$(CC) $(CFLAGS) $^ -o $@

run_read_data: read_data
	./read_data \
	../test_data/siftsmall/siftsmall_base.fvecs \
	../test_data/siftsmall/siftsmall_groundtruth.ivecs

run_func_read_data: func_read_data
	./func_read_data ../test_data/siftsmall/siftsmall_base.fvecs

clean:
	rm -f $(FILES) *.o
