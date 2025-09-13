CC?=g++
CFLAGS?=-Wall -Wextra -Wpedantic -O3
CFLAGS+=-Ilib/cglm/include
LDFLAGS?=-O3 -flto
LDFLAGS+=-lm -lSDL3

SRCS=src/main.c src/util.c src/rt.c

rt: $(SRCS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^

.PHONY: clean

clean:
	-rm rt
