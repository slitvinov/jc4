.SUFFIXES:
.SUFFIXES: .c

all: \
value \

.c:; c99 $(CFLAGS) -o $@ $< `pkg-config --libs --cflags gsl`

