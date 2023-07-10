CC = g++
CFLAGS = -Wall -std=c++14 -g
LDFLAGS = -lpthread
SRCS = main.cpp kmerpos.cpp
OBJS = $(SRCS:.cpp=.o)
TARGET = hicmap

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)
