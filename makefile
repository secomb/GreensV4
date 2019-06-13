# Place this make file in the same directory as the source files.
# Type " make clean " to remove previous build files. 
# Compile in the command line by typing:
# make
# This creates an executable "flowEstimate", which can be run by typing:
# ./flowEstimate

CC=g++ # define the compiler to use
TARGET=greens# define the name of the executable
SOURCES=analyzenet.cpp bicgstab.cpp blood.cpp cmgui.cpp contour.cpp contr_lines.cpp contr_shade.cpp convect.cpp eval.cpp gaussj.cpp greens.cpp histogram.cpp initgreens.cpp input.cpp ludcmp.cpp main.cpp nrutil.cpp outboun.cpp picturenetwork.cpp postgreens.cpp putrank.cpp setuparrays0.cpp setuparrays1.cpp setuparrays2.cpp testconvect.cpp tissrate.cpp # list source files
CFLAGS=-O3
LFLAGS=-Wall -lm -std=c++17 # -lboost_filesystem -lboost_system

# define list of objects
OBJSC=$(SOURCES:.c=.o)
OBJS=$(OBJSC:.cpp=.o)

# the target is obtained linking all .o files
# for mac: remove -lstdc++fs tag following $(TARGET)
all: $(SOURCES) $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o $(TARGET) -lstdc++fs

purge: clean
	rm -f $(TARGET)

clean:
	rm -f *.o
