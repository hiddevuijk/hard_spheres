TARGET = test.exe
OBJS = main.o pair_correlation.o xyz.o randuni.o  configFile.o
CC = g++
#CFLAGS = -c -Wall -g
#LFLAGS = -Wall -g
CFLAGS = -c -Wall -O3 -DNDEBUG -std=c++11
LFLAGS = -Wall -O3 -DNDEBUG

$(TARGET): $(OBJS)
	$(CC) $(LFLAGS)  $(OBJS) -o $(TARGET)

main.o: main.cpp xyz.h system.h pair_correlation.h configFile.h
	$(CC) $(CFLAGS) main.cpp

pair_correlation.o: pair_correlation.cpp pair_correlation.h xyz.h
	$(CC) $(CFLAGS)	pair_correlation.cpp

xyz.o: xyz.cpp xyz.h
	$(CC) $(CFLAGS) xyz.cpp

configFile.o: configFile.cpp configFile.h
	$(CC) $(CFLAGS) configFile.cpp

randuni.o: randuni.cpp randuni.h
	$(CC) $(CFLAGS) randuni.cpp


.PHONY: clean
clean:
	rm -f  $(OBJS) $(TARGET)

.PHONY: cleanObject
cleanObject:
	rm -f  $(OBJS)

