# pi.out: pi.cpp
# 	gcc -o pi.out pi.cpp

CC = g++
CFLAGS = -Wall -O2
SOURCE = pi.cpp
TARGET = pi.out

all: $(TARGET)

$(TARGET): $(SOURCE)
	$(CC) $(CFLAGS) -o $(TARGET) $(SOURCE)

clean:
	rm -f *.o *.s *.exe *.out