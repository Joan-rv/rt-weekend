CC?=g++
CFLAGS?=-Wall -Wextra -Wpedantic -O3 -ggdb
CFLAGS+=-Ilib/cglm/include -Ilib/stb
LDFLAGS+=-lm -lSDL3

SRCS=src/main.c src/util.c src/rt.c src/stb.c

rt: $(SRCS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^

.PHONY: clean

clean:
	-rm rt
