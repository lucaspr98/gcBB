CC = gcc
CFLAGS = -c -std=gnu99 
OBJFILES = main.o boss.o bwsd.o
TARGET = gcBB

# coverage information
COVERAGE = 0

DEFINES = -DCOVERAGE=$(COVERAGE)

all: $(TARGET)
	make -C egap/

$(TARGET): $(OBJFILES)
	$(CC) -g -O0 -o $(TARGET) $(OBJFILES) -lm

boss.o: boss.c boss.h
	$(CC) $(CFLAGS) boss.c -o boss.o

bwsd.o: bwsd.c bwsd.h
	$(CC) $(CFLAGS) $(DEFINES) bwsd.c -o bwsd.o

main.o: main.c
	$(CC) $(CFLAGS) main.c -o main.o

clean:
	rm -f $(TARGET) $(OBJFILES) *~
