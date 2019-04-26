OBJS = compute.o import.o matrix.o
CC = gcc
compute : $(OBJS)
	$(CC) -o compute $(OBJS)
$(OBJS) : type.h import.h
compute.o : compute.c
matrix.o : matrix.c
import.o : import.c
.PHONY : clean
clean :
	rm -fr *.o compute