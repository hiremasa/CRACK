UNAME := $(shell uname)

SRCS = main.cpp data.cpp tree.cpp forest.cpp score.cpp synthetic_data.cpp
OBJS = $(SRCS:.cpp=.o)


CC = g++
LDFLAGS=
CFLAGS= -std=c++11 -O3

# Linux addon
ifeq ($(UNAME), Linux)
LDFLAGS += -fopenmp
CFLAGS += -fopenmp
endif

all : crack.run 

crack.run: $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(OBJS)


%.o : %.cpp
	 $(CC) $(CFLAGS) -c $*.cpp

clean:
	-rm *.o crack.run

cp:
	cp crack.run ../testPackage

depend:
	makedepend -Y -- $(CFLAGS) -- $(SRCS)
# DO NOT DELETE

data.o: defs.h
score.o: tree.h score.h
tree.o: data.h tree.h defs.h
synthetic_data.o: synthetic_data.h data.h tree.h defs.h
forest.o: data.h tree.h defs.h score.h forest.h
main.o: data.h score.h forest.h tree.h rack_1.h rack_2.h Wall_Time.h synthetic_data.h
