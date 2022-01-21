CC = g++
CFLAGS = -O1 -fno-tree-vectorize -Wall -Wextra -Wpedantic # -g -pg
FILES = read_hyperplanes

.PHONY: default all clean

default: read_hyperplanes

all: $(FILES)

%: %.cpp
	$(CC) $(CFLAGS) $^ -o $@

clean:
	rm -f $(FILES) *.o
