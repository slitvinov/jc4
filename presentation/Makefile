.SUFFIXES:
.SUFFIXES: .c

all: \
value \
rl \
irl \
mapox \

CFLAGS = -O3

.c:; c99 $(CFLAGS) -o $@ $< `pkg-config --libs --cflags gsl`

