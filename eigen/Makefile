CC	= gcc
LIBS    = ./arpack/libarpack_SUN4.a -lgfortran -lm
OPTIMIZE= -O0
CFLAGS  = -c $(OPTIMIZE) -g -W -Wall -Wunused -Wcast-align -Werror -fstack-protector-all -Wfloat-equal -Wpointer-arith -Wwrite-strings -Wcast-align -Wno-long-long -Wmissing-declarations
LDFLAGS =    $(OPTIMIZE) -g -W -Wall -Wunused -Wcast-align -Werror -fstack-protector-all -Wfloat-equal -Wpointer-arith -Wwrite-strings -Wcast-align -Wno-long-long -Wmissing-declarations

HFILES	= f.h
SOURCES = main.c numsds_spectral.c s.c calc_mesh.c funcs.c
OBJECTS = $(addprefix $(OBJS_FOLDER),$(SOURCES:.c=.o))

OBJS_FOLDER = ./objs/
EXECUTABLE = a.out

all: $(HFILES) $(SOURCES) $(EXECUTABLE)
    
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ $(LIBS)

$(OBJS_FOLDER)%.o: %.c $(HFILES)
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f $(EXECUTABLE) $(OBJECTS) leak.out* log.txt
