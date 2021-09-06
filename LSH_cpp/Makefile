CC = g++
CFLAGS = -Wall -Wextra -Wpedantic
FILES = lsh

.PHONY: default all clean

default: all

all: $(FILES)

%.o: %.c
	$(CC) $(CFLAGS) -c $^

clean:
	rm -f $(FILES) *.o
