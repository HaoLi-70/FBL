# Makefile MCMC inversion
#-------------------------------------------
CC = gcc
ARRCS = ar rcs
CFLAGS = -c


all: Run libmine.a clean

OBJ1=ALLOCATION.o LL04.o LU.o RDATOM.o RDINPUT.o\
        RDMODEL.o STRING2.o TIME_PRINT.o\
       
OBJ2=main.o	

ALLOCATION.o:ALLOCATION.c
	$(CC) $(CFLAGS) ALLOCATION.c

LL04.o:LL04.c 
	$(CC) $(CFLAGS) LL04.c

LU.o:LL04.c 
	$(CC) $(CFLAGS) LU.c

RDATOM.o:RDATOM.c
	$(CC) $(CFLAGS) RDATOM.c

RDINPUT.o:RDATOM.c
	$(CC) $(CFLAGS) RDINPUT.c

RDMODEL.o:RDATOM.c
	$(CC) $(CFLAGS) RDMODEL.c

STRING2.o:STRING2.c
	$(CC) $(CFLAGS) STRING2.c

TIME_PRINT.o:STRING2.c
	$(CC) $(CFLAGS) TIME_PRINT.c

libmine.a:$(OBJ1)
	$(ARRCS) libmine.a $(OBJ1)

main.o:main.c
	$(CC) $(CFLAGS) main.c

Run: $(OBJ2) libmine.a
	$(CC) $(OBJ2) -o Run -L. -lm -lmine
	mv Run ../

.PHONY.:clean

clean:
	rm -f *.o *.a


