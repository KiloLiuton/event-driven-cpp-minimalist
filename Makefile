# Usage: make N=100 K=10 p=0.1
# Creates a topology file in the folder TOPOLOGIES_DIR which will be used to compile an executable
# for simulations of the WCM dynamics under the generated topology. Executable files are put into
# the bin folder.

SHELL := /bin/bash
ODIR = obj
EXEC_DIR = bin
SRC_DIR = src
INCLUDE_DIR = include
TOPOLOGIES_DIR = topology_headers

_OBJ = main.o dynamics.o
OBJ = $(patsubst %, $(ODIR)/%, $(_OBJ))

N ?= 20
K ?= 3
p ?= 0.000000
topology_seed ?= 42
N_ = $(shell printf %05d $(N))
K_ = $(shell printf %04d $(K))
p_ = $(shell echo $(p) | sed 's/\./_/')
s_ = $(shell printf %d $(topology_seed))
topology_file = "$(TOPOLOGIES_DIR)/$(N_)-$(K_)-$(p_)-gseed_$(s_).hpp"
prog = "$(EXEC_DIR)/sim-$(N_)-$(K_)-$(p_)-gseed_$(s_)"
# prog = "sim-" + $(topology_file)

CC = g++
#CFLAGS = -Wall -std=c++11 -O3 -march=native -fopenmp -include $(topology_file) -I$(INCLUDE_DIR) -I$(TOPOLOGIES_DIR) -Ipcg_random/ -Iexperiments
CFLAGS = -g -Wall -std=c++11 -fopenmp -include $(topology_file) -I$(INCLUDE_DIR) -I$(TOPOLOGIES_DIR) -Ipcg_random/ -Iexperiments
LDFLAGS = -lm -fopenmp

# all: $(prog) compile_commands.json
all: $(prog)

$(topology_file):
	test -x $(EXEC_DIR)/build_header || { $(MAKE) compile_header_builder;}
	./$(EXEC_DIR)/build_header $(N) $(K) $(p) $(TOPOLOGIES_DIR) $(topology_seed)

# For each .o file. '$@' = left of :, '$<' = first dependency (the .cpp file)
$(ODIR)/%.o: $(SRC_DIR)/%.cpp $(topology_file)
	$(CC) -c -o $@ $< $(CFLAGS)

# Link all objects. '$^' = all dependencies
$(prog): $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

# $compile_commands.json: Makefile
# 	bear $(MAKE)

compile_header_builder:
	$(CC) $(SRC_DIR)/build_header.cpp -Wall -std=c++11 -O3 -march=native -o $(EXEC_DIR)/build_header

.PHONY: clean
clean:
	rm -f $(EXEC_DIR)/sim-* $(ODIR)/*.o
