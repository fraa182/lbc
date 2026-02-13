CC = gcc
CFLAGS = -Wall -O3 -I../../include
LDFLAGS = -lm -fopenmp

# Source directories
SRC_DIR = ../../src
CORE_DIR = $(SRC_DIR)/core
BC_DIR = $(SRC_DIR)/bc
UTILS_DIR = $(SRC_DIR)/utils

# Source files
SRCS = $(wildcard *.c) \
       $(wildcard $(CORE_DIR)/*.c) \
       $(wildcard $(BC_DIR)/*.c) \
       $(wildcard $(UTILS_DIR)/*.c)

# Target executable
TARGET = couette

.PHONY: all clean

all: $(TARGET).o

$(TARGET).o: $(SRCS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm -f $(TARGET).o
	rm -f sol/*.vtk
