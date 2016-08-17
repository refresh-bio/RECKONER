KMC_DIR = kmc_dir
RECKONER_DIR = reckoner_dir
BIN_DIR = bin

export KMC_DIR

all: kmc_main reckoner

kmc_main:
	cd $(KMC_DIR) && $(MAKE) DISABLE_ASMLIB=true
	-cp $(KMC_DIR)/bin/kmc $(BIN_DIR)
	-cp $(KMC_DIR)/bin/kmc_tools $(BIN_DIR)

reckoner:
	cd $(RECKONER_DIR) && $(MAKE)
	-cp $(RECKONER_DIR)/reckoner $(BIN_DIR)

clean:
	-rm -f $(BIN_DIR)/cutter
	cd $(CUTTER_DIR) && $(MAKE) clean
	-rm -f $(BIN_DIR)/reckoner
	cd $(RECKONER_DIR) && $(MAKE) clean
	-rm -f $(BIN_DIR)/kmc_tools
	-rm -f $(BIN_DIR)/kmc
	cd $(KMC_DIR) && $(MAKE) clean
