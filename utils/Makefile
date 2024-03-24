CC = gcc
CFLAGS = -g -Wall -Wno-char-subscripts -Wno-unused-function -std=gnu11
LIBS = -lm #-ldl 
MCOBJ =
USEMC = 0

SRC = $(wildcard *.c)

nj: $(SRC:%.c=%.o)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

clean:
	\rm *.o nj


