CC = g++
CFLAGS = -std=c++11

SRC_DIR = src
INCLUDE_DIR = include

# List of source files
SRCS = $(SRC_DIR)/main.cpp $(SRC_DIR)/read.cpp $(SRC_DIR)/solve.cpp $(SRC_DIR)/objectfunc.cpp $(SRC_DIR)/initialize.cpp $(SRC_DIR)/display.cpp $(SRC_DIR)/constructor.cpp $(SRC_DIR)/constraints.cpp $(SRC_DIR)/node.cpp $(SRC_DIR)/edge.cpp $(SRC_DIR)/element.cpp

# List of header files
HEADERS = $(INCLUDE_DIR)/read.h $(INCLUDE_DIR)/solve.h $(INCLUDE_DIR)/objectfunc.h $(INCLUDE_DIR)/initialize.h $(INCLUDE_DIR)/display.h $(INCLUDE_DIR)/constructor.h $(INCLUDE_DIR)/constraints.h $(INCLUDE_DIR)/node.h $(INCLUDE_DIR)/edge.h $(INCLUDE_DIR)/element.h

# Output executable
TARGET = myprogram

# Default target
all: $(TARGET)

# Rule to build the executable
$(TARGET): $(SRCS) $(HEADERS)
	$(CC) $(CFLAGS) -I$(INCLUDE_DIR) -o $(TARGET) $(SRCS) -lnlopt

# Clean rule
clean:
	rm -f $(TARGET)

.PHONY: all clean
