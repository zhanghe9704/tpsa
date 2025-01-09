CC = g++
INSTALL_DIR = ..
CFLAGS = -O3 -Wall -std=c++14 

TPSA_LIBS = tpsa
LIBS = -L$(INSTALL_DIR)/lib -l$(TPSA_LIBS) -lstdc++ -lm
INCLUDES = -I$(INSTALL_DIR)/include/



SRC = $(wildcard *.cc)
TARGET = $(SRC:.cc=)

.PHONY: all
all: $(TARGET)

.cc:
	$(CC) -o $@ $< $(INCLUDES) $(LIBS)  $(CFLAGS) 

.PHONY: clean
clean:
	rm -f $(TARGET)