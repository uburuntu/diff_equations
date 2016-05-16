CFLAGS=-c -W -Wall -Wunused -Wcast-align -Werror -fstack-protector-all -Wfloat-equal -Wpointer-arith -Wwrite-strings -Wcast-align -Wno-format -Wno-long-long -Wmissing-declarations
LDFLAGS=-W -Wall -Wunused -Wcast-align -Werror -fstack-protector-all -Wfloat-equal -Wpointer-arith -Wwrite-strings -Wcast-align -Wno-format -Wno-long-long -Wmissing-declarations -lm
# You can add -Ofast flag for optimization, but che-to ne to schitaet

SOURCES=func.c gas_two.c setka.c gnuplot.c shema.c tabtex.c
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=a.out
REPORT_FILENAME=theplot.pdf

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	gcc $(OBJECTS) ./laspack/liblaspack.a -o $@ $(LDFLAGS)

.cpp.o:
	gcc $(CFLAGS) $<

clean:
	rm -rf *.o a.out leak.out output.txt

gnuplot.o: gnuplot.h

func.o: func.h

tabtex.o: tabtex.h

pdf:
	./a.out
	pdflatex -shell-escape theplot.tex
	cp ./$(REPORT_FILENAME) ./report/ --force
	mv -t ./report/report_files *.tex *.eps *.log *.gnuplot *.pdf *.aux
	cp ./report/$(REPORT_FILENAME) ./

cleanall:
	rm -rf *.o a.out leak.out output.txt *.tex *.log theplot*.*

cleanpdf:
	rm -rf *.tex *.log theplot*.*
