LEAKS = -fsanitize=address -fsanitize=leak -fsanitize=undefined -fsanitize=null -fsanitize=bounds-strict
CPPFLAGS = -O3 -g -lpthread --param large-function-growth=10000 -finline-limit=10000 --param inline-unit-growth=100000 -Winline -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format

all: a.out
a.out: $(patsubst %.cpp, %.o,$(wildcard *.cpp))
	g++ $^ $(CPPFLAGS) -o $@
%.o: %.cpp *.h
	g++ -c $< $(CPPFLAGS)
report: report.tex #*.png
	pdflatex report.tex
clean:
	rm -f *.out *.o *.zip
ready: clean
	zip Borisenkov_NN.zip *.h *.cpp Makefile
