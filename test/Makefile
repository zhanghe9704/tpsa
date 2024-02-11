#CC = g++
CC = clang++
INCDIR = ../include
CFLAGS = -O3 -Wall -std=c++14 

TPSA_LIBS = tpsa
LIBS = -L.. -l$(TPSA_LIBS) -lstdc++ -lm

SRC = $(wildcard *.cc)
TARGET = $(SRC:.cc=)
TARGET := $(filter-out catch_main, $(TARGET))
OBJ = $(SRC:.cc=.o)
OBJ2 = $(filter-out catch_main.o, $(OBJ))
DEPS =$(INCDIR)/da.h catch.hpp

#$(info    OBJ is $(OBJ))
#$(info    OBJ2 is $(OBJ2))
#$(info    TARGET is $(TARGET))


.PHONY: all
all: $(TARGET) 

$(TARGET): $(OBJ)
	$(CC) -o $@ $^ $(LIBS)  $(CFLAGS)
	
%.o: %.cc
	$(CC) -c $< 

#$(CC) -o $@ $< $(LIBS)  $(CFLAGS)		

.PHONY: clean
clean:
	rm -f $(TARGET) $(OBJ2)