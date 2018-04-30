CC = gcc

CFLAGS = -Wall -mmmx -msse -msse2 -mssse3 -msse4 -mavx2 

all: lab2.o

lab2.o: lab2.c 
	$(CC) $(CFLAGS) lab2.c -o lab2 -lm

clean:
	rm -f *.o 

