CFLAGS  = -Wall -Wextra -Werror -pedantic
SOURCES = $(wildcard *.c)
TARGETS = $(patsubst %.c,%,$(SOURCES))

all: $(TARGETS)

$(TARGETS): % : %.o
	$(CC) $(CFLAGS) -o $@ $^

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *.o $(TARGETS)

.PHONY: all clean
