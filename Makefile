CUTTER_DIR = cutter_dir
DB_CUTTER_DIR = kmc_tools_dir/kmc1_db_cutter
KMC_CONVERTER_DIR = kmc_tools_dir/kmc2DBkmc1DB_converter
KMC_DIR = kmc_dir
RECKONER_DIR = reckoner_dir
BIN_DIR = bin

export KMC_DIR

all: kmc_main reckoner cutter

kmc_main:
	cd $(KMC_DIR) && $(MAKE) DISABLE_ASMLIB=true
	-cp $(KMC_DIR)/bin/kmc $(BIN_DIR)
	-cp $(KMC_DIR)/bin/kmc_tools $(BIN_DIR)

reckoner:
	cd $(RECKONER_DIR) && $(MAKE)
	-cp $(RECKONER_DIR)/reckoner $(BIN_DIR)

cutter:
	cd $(CUTTER_DIR) && $(MAKE)
	-cp $(CUTTER_DIR)/cutter $(BIN_DIR)

clean:
	-rm -f $(BIN_DIR)/cutter
	cd $(CUTTER_DIR) && $(MAKE) clean
	-rm -f $(BIN_DIR)/reckoner
	cd $(RECKONER_DIR) && $(MAKE) clean
	-rm -f $(BIN_DIR)/kmc_tools
	-rm -f $(BIN_DIR)/kmc
	cd $(KMC_DIR) && $(MAKE) clean
