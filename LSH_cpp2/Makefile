CC = g++
CFLAGS = -Wall -Wextra -Wpedantic
FILES = lsh mlsh

.PHONY: default all clean

default: lsh

all: $(FILES)

%.o: %.c
	$(CC) $(CFLAGS) -c $^

clean:
	rm -f $(FILES) *.o
