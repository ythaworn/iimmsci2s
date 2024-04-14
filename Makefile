CC = gcc
CPPFLAGS = -O3
HDRS = paml.h

EXE = iimmsci2s
SRCS = iimmsci2s.c tools.c

LIBS = -lm

OBJS = $(SRCS:.c=.o)

all: $(EXE)

$(EXE): $(OBJS) $(HDRS)
	$(CC) $(CPPFLAGS) -o $@ $(OBJS) $(LIBS)

$(OBJS): $(HDRS)

clean:
	rm -f core $(EXE) *.o
