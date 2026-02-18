CXX = g++
CXXFLAGS = -std=c++14 -static-libstdc++ -Wall -O3 -fopenmp -I./src -w
LDFLAGS = -L./lib -lz -lspoa
LIBDIR := -L.

SRCS1 = ./src/paf.cpp ./src/mahit.cpp ./src/aln2tr.cpp ./src/trcss.cpp ./src/edlib.cpp ./src/main.cpp
OBJS1 = $(SRCS1:.cpp=.o)
TARGET1 = r2rtr

all: $(TARGET1) 

$(TARGET1): $(OBJS1)
	$(CXX) $(CXXFLAGS) $(OBJS1) -o $(TARGET1) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS1) $(OBJS2) $(TARGET1) $(TARGET2)