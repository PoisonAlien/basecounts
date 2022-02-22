CC ?= gcc
CFLAGS ?= -Wall -g -O3 -pthread

basecounts: src/basecounts.c
		$(CC) $(CFLAGS) src/basecounts.c -lhts -o basecounts
