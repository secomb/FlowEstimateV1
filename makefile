# Place this make file in the same directory as the source files.
# Type " make clean " to remove previous build files. 
# Compile in the command line by typing:
# make
# This creates an executable "flowEstimate", which can be run by typing:
# ./flowEstimate

CC=g++ # define the compiler to use
TARGET=flowEstimate # define the name of the executable
SOURCES=Amatrix.cpp Amultiply.cpp analyzenet.cpp analyzeresults.cpp bvector.cpp cmgui.cpp dishem.cpp flow.cpp flowdirections.cpp histogram.cpp input.cpp ludcmp.cpp main.cpp nrutil.cpp picturenetwork.cpp putrank.cpp setuparrays1.cpp sparse_cgsymm.cpp viscor.cpp writeflow.cpp # list source files
CFLAGS=-O3
LFLAGS=-Wall -lm 

# define list of objects
OBJSC=$(SOURCES:.c=.o)
OBJS=$(OBJSC:.cpp=.o)

# the target is obtained linking all .o files
all: $(SOURCES) $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o $(TARGET) -lstdc++fs

purge: clean
	rm -f $(TARGET)

clean:
	rm -f *.o
