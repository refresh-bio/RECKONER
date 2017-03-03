KMC_DIR = kmc_dir/
RECKONER_DIR = reckoner_dir/
BIN_DIR = bin/
CC=g++

all: kmc_main reckoner

kmc_main:
	cd $(KMC_DIR) && $(MAKE) DISABLE_ASMLIB=true CC=$(CC)
	-cp $(KMC_DIR)bin/kmc $(BIN_DIR)
	-cp $(KMC_DIR)bin/kmc_tools $(BIN_DIR)

reckoner: kmc_main
	cd $(RECKONER_DIR) && $(MAKE) KMC_DIR=$(KMC_DIR) CC=$(CC)
	-cp $(RECKONER_DIR)reckoner $(BIN_DIR)

clean:
	-rm -f $(BIN_DIR)*
	cd $(RECKONER_DIR) && $(MAKE) clean
	cd $(KMC_DIR) && $(MAKE) clean
