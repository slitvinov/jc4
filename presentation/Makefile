.SUFFIXES:
.SUFFIXES: .c

all: \
value \
rl \

.c:; c99 $(CFLAGS) -o $@ $< `pkg-config --libs --cflags gsl`

