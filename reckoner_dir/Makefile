KMC_DIR_IN = ../$(KMC_DIR)

CC = g++
CFLAGS = -Wall -I$(KMC_DIR_IN) -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archiv -Ilibs -O3 -fopenmp -std=c++11
LIBS = ../kmc_dir/3rd_party/cloudflare/libz.a
SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:.cpp=.o)

all: reckoner

reckoner: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(KMC_DIR_IN)/kmc_api/*.o $(LIBS) -o $@

$(OBJS): %.o : %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	-rm -f reckoner $(OBJS)
