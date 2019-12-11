TARGET = test.exe
OBJS = main.o pair_correlation.o xyz.o system.o configFile.o neighbour_check.o
CC = g++
#CFLAGS = -c -Wall -g -std=c++11
#LFLAGS = -Wall -g
CFLAGS = -c -Wall -O3 -DNDEBUG -std=c++11
LFLAGS = -Wall -O3 -DNDEBUG

$(TARGET): $(OBJS)
	$(CC) $(LFLAGS)  $(OBJS) -o $(TARGET)

pair_correlation.o: pair_correlation.cpp pair_correlation.h xyz.h
	$(CC) $(CFLAGS)	pair_correlation.cpp

xyz.o: xyz.cpp xyz.h
	$(CC) $(CFLAGS) xyz.cpp

system.o: system.cpp system.h 
	$(CC) $(CFLAGS) system.cpp

configFile.o: configFile.cpp configFile.h
	$(CC) $(CFLAGS) configFile.cpp

main.o: main.cpp xyz.h system.h pair_correlation.h  potential.h neighbour_check.h
	$(CC) $(CFLAGS) main.cpp

neighbour_check.o: neighbour_check.h neighbour_check.cpp
	$(CC) $(CFLAGS) neighbour_check.cpp


.PHONY: clean
clean:
	rm -f  $(OBJS) $(TARGET) 

.PHONY: cleanObject
cleanObject:
	rm -f  $(OBJS)

