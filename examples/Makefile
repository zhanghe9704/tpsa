#CC = g++
CC = clang++
INCDIR = ../include
CFLAGS = -O3 -Wall -std=c++14 

TPSA_LIBS = tpsa
LIBS = -L.. -l$(TPSA_LIBS) -lstdc++ -lm

SRC = $(wildcard *.cc)
TARGET = $(SRC:.cc=)
DEPS =$(INCDIR)/da.h

.PHONY: all
all: $(TARGET)

.cc:
	$(CC) -o $@ $< $(LIBS)  $(CFLAGS) 

.PHONY: clean
clean:
	rm -f $(TARGET)