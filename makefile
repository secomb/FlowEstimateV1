# Updated makefile for flow estimation algorithm by Fryet al.
# Compile in command line by executing 'make', which creates executable "flowEstimate".
# Last updated by Grace Lee, September 2020

# Compiler options
CXX = g++
TARGET = flowEstimate

# Specify object build directory and sources
OBJ_DIR = ./obj
SRC = $(wildcard *.cpp)
OBJS = $(SRC:%.cpp=$(OBJ_DIR)/%.o)

# Compiler flags
CXXFLAGS = -std=c++17
LDFLAGS =

all: build $(TARGET)

# Compile cpp files to objects in ./obj
$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(LDFLAGS)

# Build target by compiling .o files
$(TARGET): $(OBJS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $^ $(LDFLAGS)

.PHONY: all build clean debug release

build:
	@mkdir -p $(OBJ_DIR)

debug: CXXFLAGS += -g
debug: all

release: CXXFLAGS += -O3
release: all

clean:
	-@rm -rvf $(OBJ_DIR)/*
	-@rm -vf $(TARGET)
