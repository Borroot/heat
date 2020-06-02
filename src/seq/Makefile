CFLAGS  = -Wall -Wextra -Werror -pedantic -O2 -fopenmp
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
