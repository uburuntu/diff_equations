CFLAGS=-c -Ofast -W -Wall -Wunused -Wcast-align -Werror -fstack-protector-all -Wfloat-equal -Wpointer-arith -Wwrite-strings -Wcast-align -Wno-format -Wno-long-long -Wmissing-declarations
LDFLAGS=-Ofast -W -Wall -Wunused -Wcast-align -Werror -fstack-protector-all -Wfloat-equal -Wpointer-arith -Wwrite-strings -Wcast-align -Wno-format -Wno-long-long -Wmissing-declarations -lm
SOURCES=func.c gas_two.c setka.c shema.c tabtex.c
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=a.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	gcc $(OBJECTS) ./laspack/liblaspack.a -o $@ $(LDFLAGS)

.cpp.o:
	gcc $(CFLAGS) $<

clean:
	rm -rf *.o a.out leak.out output.txt

func.o: func.h

tabtex.o: tabtex.h