#============================================================================================
# definitions
#============================================================================================
TOOLS_DIR = ./Tools

MAKE_COMMAND = make -f Makefile.in TOOLS_DIR=$(TOOLS_DIR)

#============================================================================================
# rules
#============================================================================================
all: 
	$(MAKE_COMMAND) all

lib: 
	$(MAKE_COMMAND) lib

bin: 
	$(MAKE_COMMAND) bin

#============================================================================================
# rules to clean up
#============================================================================================
clean:
	rm -f *.o
	rm -f ./CSpack-Tools/*.o

cleanall:
	rm -f *.o *~ CSpack_main libCSpack.a
	rm -f ./CSpack-Tools/*.o ./CSpack-Tools/*.~

cleanallDEV:
	make cleanall
	rm -f $(TOOLS_DIR)/*.o $(TOOLS_DIR)/*.~
	rm -f $(TOOLS_DIR)/Definitions/*.o $(TOOLS_DIR)/Definitions/*.~
	rm -f $(TOOLS_DIR)/Simple_routines/*.o $(TOOLS_DIR)/Simple_routines/*.~
	rm -f $(TOOLS_DIR)/Cosmology/*.o $(TOOLS_DIR)/Cosmology/*.~
	rm -f $(TOOLS_DIR)/Integration/*.o $(TOOLS_DIR)/Integration/*~

tidy:
	make cleanallDEV

wipeDS:
	find . -type f -name \.DS_Store -print | xargs rm

#============================================================================================
