# pi.out: pi.cpp
# 	gcc -o pi.out pi.cpp

CC = g++
CFLAGS = -Wall -O2
LDFLAGS = -pthread
SOURCE = pi.cpp
TARGET = pi.out

all: $(TARGET)

$(TARGET): $(SOURCE)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $(TARGET) $(SOURCE)

clean:
	rm -f *.o *.s *.exe *.out