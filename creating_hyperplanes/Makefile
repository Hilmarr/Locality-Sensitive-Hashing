CC = g++
CFLAGS = -O1 -fno-tree-vectorize -Wall -Wextra -Wpedantic # -g -pg
FILES = read_hyperplanes creating_hplanes

.PHONY: default all clean

default: all

all: $(FILES)

%: %.cpp
	$(CC) $(CFLAGS) $^ -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

creating_hplanes: creating_hplanes.o read_vector_file_funcs.o
	$(CC) $(CFLAGS) $^ -o $@

run_creating_hplanes: creating_hplanes
	./creating_hplanes \
	../test_data/siftsmall/siftsmall_learn.fvecs

clean:
	rm -f $(FILES) *.o
