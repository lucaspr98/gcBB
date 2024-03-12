CC = gcc
CFLAGS = -O3 -Wall -Wno-char-subscripts -Wno-unused-function -c -std=gnu99 
#CFLAGS = -g -O0
OBJFILES = external.o boss.o bwsd.o lib/rankbv.o
TARGET = gcBB

COVERAGE = 0
ALL_VS_ALL = 1
DEBUG = 0

DEFINES = -DCOVERAGE=$(COVERAGE) -DALL_VS_ALL=$(ALL_VS_ALL) -DDEBUG=$(DEBUG)

all: $(TARGET)
	make -C egap/ && make -C utils/

$(TARGET): main.c $(OBJFILES) 
	$(CC) $^ -o $(TARGET) -ldl -lm 

%.o: %.c %.h
	$(CC) $(CFLAGS) $(DEFINES) -c $< -o $@

clean:
	rm -f $(TARGET) $(OBJFILES) *~ && cd utils && rm *.o 
