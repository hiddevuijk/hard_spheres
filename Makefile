TARGET = test.exe
OBJS = main.o pair_correlation.o xyz.o system.o
CC = g++
#lCFLAGS = -c -Wall -g -std=c++11
#lLFLAGS = -Wall -g
CFLAGS = -c -Wall -O3 -DNDEBUG -std=c++11
LFLAGS = -Wall -O3 -DNDEBUG

$(TARGET): $(OBJS)
	$(CC) $(LFLAGS)  $(OBJS) -o $(TARGET)

main.o: main.cpp xyz.h system.h pair_correlation.h  potential.h
	$(CC) $(CFLAGS) main.cpp

pair_correlation.o: pair_correlation.cpp pair_correlation.h xyz.h
	$(CC) $(CFLAGS)	pair_correlation.cpp

xyz.o: xyz.cpp xyz.h
	$(CC) $(CFLAGS) xyz.cpp

system.o: system.cpp system.h 
	$(CC) $(CFLAGS) system.cpp



.PHONY: clean
clean:
	rm -f  $(OBJS) $(TARGET)

.PHONY: cleanObject
cleanObject:
	rm -f  $(OBJS)

