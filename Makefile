CC?=g++
CFLAGS?=-Wall -Wextra -Wpedantic -O3
CFLAGS+=-Ilib/cglm/include
LD=$(CC)
LDFLAGS?=-O3 -flto
LDFLAGS+=-lm -lSDL3

OBJS=main.o

rt: $(OBJS)
	$(LD) $(LDFLAGS) -o $@ $^

.PHONY: clean fclean
clean:
	-rm -f $(OBJS)

fclean: clean
	-rm rt
