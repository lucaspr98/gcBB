CC = gcc
CFLAGS = -O3 -Wall -Wno-char-subscripts -Wno-unused-function -c -std=gnu99 
OBJFILES = main.o boss.o bwsd.o lib/rankbv.o
TARGET = gcBB

COVERAGE = 0
BWSD_ALL = 1

DEFINES = -DCOVERAGE=$(COVERAGE) -DBWSD_ALL=$(BWSD_ALL)

all: $(TARGET)
	make -C egap/ && make -C utils/

$(TARGET): $(OBJFILES) 
	$(CC) -o $(TARGET) $(OBJFILES) -ldl -lm 

boss.o: boss.c boss.h
	$(CC) $(CFLAGS) $(DEFINES) boss.c -o boss.o

bwsd.o: bwsd.c bwsd.h
	$(CC) $(CFLAGS) $(DEFINES) bwsd.c -o bwsd.o

main.o: main.c
	$(CC) $(CFLAGS) $(DEFINES) main.c -o main.o

clean:
	rm -f $(TARGET) $(OBJFILES) *~ && cd utils && rm *.o 