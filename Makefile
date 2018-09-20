ODIR = obj

_OBJ = main.o dynamics.o
OBJ = $(patsubst %, $(ODIR)/%, $(_OBJ))

CC = g++
CFLAGS = -Wall -std=c++11 -O3 -march=native -I pcg_random -fopenmp
LDFLAGS = -lm -fopenmp

N ?= 20
K ?= 3
p ?= 0.000000
s ?= 42
N_ = $(shell printf %05d $(N))
K_ = $(shell printf %04d $(K))
p_ = $(shell echo $(p) | sed 's/\./_/')
s_ = $(shell printf %d $(s))
h = "$(N_)-$(K_)-$(p_)-seed_$(s_).hpp"
prog = "sim-$(N_)-$(K_)-$(p_)-graphseed_$(s_)"

all: $(prog) compile_commands.json

$(h):
	./build_header.py $(N) $(K) $(p) $(s)

# For each .o file. '$@' = left of :, '$<' = first dependency (the .c file)
$(ODIR)/%.o: %.cpp $(h)
	$(CC) -include $(h) -c -o $@ $< $(CFLAGS)

# Link all objects. '$^' = all dependencies
$(prog): $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

$compile_commands.json: Makefile
	bear $(MAKE)

.PHONY: clean
clean:
	rm -f $(prog) $(ODIR)/*.o
