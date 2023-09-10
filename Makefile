CC = gcc
CFLAGS = -O -Wall -Wno-char-subscripts -Wno-unused-function -c -std=gnu99 
OBJFILES = main.o  external.o boss.o bwsd.o lib/rankbv.o
TARGET = gcBB

COVERAGE = 0
ALL_VS_ALL = 0
DEBUG = 0

DEFINES = -DCOVERAGE=$(COVERAGE) -DALL_VS_ALL=$(ALL_VS_ALL) -DDEBUG=$(DEBUG)

all: $(TARGET)
	make -C egap/ && make -C utils/

$(TARGET): $(OBJFILES) 
	$(CC) -o $(TARGET) $(OBJFILES) -ldl -lm 

external.o: external.c external.h
	$(CC) $(CFLAGS) $(DEFINES) external.c -o external.o

boss.o: boss.c boss.h
	$(CC) $(CFLAGS) $(DEFINES) boss.c -o boss.o

bwsd.o: bwsd.c bwsd.h
	$(CC) $(CFLAGS) $(DEFINES) bwsd.c -o bwsd.o

main.o: main.c
	$(CC) $(CFLAGS) $(DEFINES) main.c -o main.o

clean:
	rm -f $(TARGET) $(OBJFILES) *~ && cd utils && rm *.o 