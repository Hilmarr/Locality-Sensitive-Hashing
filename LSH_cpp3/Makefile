CC = pgc++
CFLAGS = -Wall -Wextra -Wpedantic -acc -fast -Minfo=all
FILES = lsh

.PHONY: default all clean profile_lsh

default: lsh

all: $(FILES)

%: %.cpp
	$(CC) $(CFLAGS) $^ -o $@

profile_lsh:
	gprof lsh gmon.out > gprof_results.txt

nvprof_lsh:
	nvprof --unified-memory-profiling off ./lsh


clean:
	rm -f $(FILES) *.o
