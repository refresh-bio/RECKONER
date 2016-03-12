KMC_DIR_IN = ../$(KMC_DIR)

CC=g++
CFLAGS = -Wall -I$(KMC_DIR_IN) -Ilibs -O3 -fopenmp -std=c++11
LIBS = libs/libz.a
SRCS=$(wildcard *.cpp)
OBJS=$(SRCS:.cpp=.o)

all: reckoner kmc_api

kmc_api:
	cd $(KMC_DIR_IN) && $(MAKE) kmc_api

reckoner: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(KMC_DIR_IN)/kmc_api/*.o $(LIBS) -o $@

$(OBJS): %.o : %.cpp
	$(CC) $(CFLAGS) $(OPT) -c $< -o $@

clean:
	-rm -f reckoner $(OBJS)
