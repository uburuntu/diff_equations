CC	= gcc
HFILES	= ./f.h  
OBJS	= ./objs/main.o ./objs/numsds_spectral.o ./objs/s.o ./objs/calc_mesh.o ./objs/funcs.o
LIBS    = -lgfortran
G	= -O  -g -Wall 


./m.e:	$(OBJS)
	$(CC) $(OBJS) ./arpack/libarpack_SUN4.a -lm $(LIBS)  -o ./m.e $(G)$ 

./objs/main.o:	./main.c $(HFILES)
	$(CC) $(G)$ -c ./main.c -o ./objs/main.o
./objs/numsds_spectral.o:	./numsds_spectral.c $(HFILES)
	$(CC) $(G)$ -c ./numsds_spectral.c -o ./objs/numsds_spectral.o
./objs/s.o:	./s.c $(HFILES)
	$(CC) $(G)$ -c ./s.c -o ./objs/s.o
./objs/calc_mesh.o:	./s.c $(HFILES)
	$(CC) $(G)$ -c ./calc_mesh.c -o ./objs/calc_mesh.o
./objs/funcs.o:	./s.c $(HFILES)
	$(CC) $(G)$ -c ./funcs.c -o ./objs/funcs.o

clean:
	rm -f ./m.e ./objs/*.o ./numres/* ./*.dat


