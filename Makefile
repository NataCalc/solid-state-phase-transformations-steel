CC = g++
CFLAGS = -g -O3 -I/opt/local/include -pthread
COMPILE = $(CC) $(CFLAGS) -c
OBJFILES := $(patsubst %.cpp,%.o,$(wildcard *.cpp))

all: preci

preci: $(OBJFILES)
	$(CC) -o preci $(OBJFILES)

%.o: %.cpp
	$(COMPILE) -o $@ $<

clean:
	rm *.o && rm preci