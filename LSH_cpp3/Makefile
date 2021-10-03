# CC = pgc++
# # CFLAGS = -Wall -Wextra -Wpedantic -acc -fast -Minfo=all
CC = g++
# CFLAGS = -Wall -Wextra -Wpedantic
# CFLAGS = -O3 -Wall -Wextra -Wpedantic
CFLAGS = -O1 -fno-tree-vectorize -Wall -Wextra -Wpedantic -g -pg
FILES = lsh

.PHONY: default all clean profile_lsh

default: lsh

all: $(FILES)

%: %.cpp
	$(CC) $(CFLAGS) $^ -o $@

profile_lsh:
	gprof lsh gmon.out > gprof_results.txt

clean:
	rm -f $(FILES) *.o
