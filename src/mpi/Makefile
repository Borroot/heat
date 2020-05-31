CFLAGS  = -Wall -Wextra -pedantic -O3
CC      = mpicc

SOURCES = $(shell find . -type f -name '*.c')
OBJECTS = $(patsubst %.c, %.o, $(SOURCES))
TARGETS = $(patsubst %.c, %,   $(SOURCES))

all: $(TARGETS)

$(TARGETS): % : %.o
	$(CC) $(CFLAGS) -o $@ $^

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJECTS) $(TARGETS)

.PHONY: all clean
