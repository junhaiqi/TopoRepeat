CXX = g++
CXXFLAGS = -std=c++14 -static-libstdc++ -Wall -O3 -fopenmp -I. -w
LDFLAGS = -L. -lz
LIBDIR := -L.

SRCS1 = edlib.cpp clust.cpp main.cpp
OBJS1 = $(SRCS1:.cpp=.o)
TARGET1 = raEDClust

all: $(TARGET1) 

$(TARGET1): $(OBJS1)
	$(CXX) $(CXXFLAGS) $(OBJS1) -o $(TARGET1) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS1) $(OBJS2) $(TARGET1) $(TARGET2)