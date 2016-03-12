KMC_DIR_IN = ../$(KMC_DIR)

CC = g++
CFLAGS = -Wall -O3 -std=c++11 -I$(KMC_DIR_IN)
LDFLAGS = -L$(KMC_DIR_IN)
SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:.cpp=.o)

all: kmc_api cutter

kmc_api:
	cd $(KMC_DIR_IN) && $(MAKE) kmc_api

cutter: $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(OBJS) $(KMC_DIR_IN)/kmc_api/*.o -o $@

$(OBJS): %.o : %.cpp
	$(CC) $(CFLAGS) -I$(KMC_DIR_IN) -c $< -o $@

clean:
	cd $(KMC_DIR_IN) && $(MAKE) clean
	-rm -f cutter *.o

