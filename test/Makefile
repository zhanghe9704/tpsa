CC = g++

INSTALL_DIR = ..
#CFLAGS = -O3 -Wall -std=c++14 
CFLAGS = -g -Wall -std=c++14 

TPSA_LIBS = tpsa
LIBS = -L$(INSTALL_DIR)/lib -l$(TPSA_LIBS) -lgcc -lstdc++ -lm

INCLUDES = -I$(INSTALL_DIR)/include/

SRC = $(wildcard *.cc)
TARGET = $(SRC:.cc=)
TARGET := $(filter-out catch_main, $(TARGET))
OBJ = $(SRC:.cc=.o)
OBJ2 = $(filter-out catch_main.o, $(OBJ))
#DEPS =$(INCDIR)/da.h catch.hpp

#$(info    OBJ is $(OBJ))
#$(info    OBJ2 is $(OBJ2))
#$(info    TARGET is $(TARGET))


.PHONY: all
all: $(TARGET) 

$(TARGET): $(OBJ)
	$(CC) -o $@ $^ $(INCLUDES)  $(CFLAGS)  $(LIBS)
	
%.o: %.cc
	$(CC) -c $< 	

.PHONY: clean
clean:
	rm -f $(TARGET) $(OBJ2)
.PHONY: clean_all
clean_all:
	rm -f $(TARGET) $(OBJ)