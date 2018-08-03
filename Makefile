ODIR = obj

_OBJ = main.o
OBJ = $(patsubst %, $(ODIR)/%, $(_OBJ))

CC = g++
CFLAGS = -Wall -std=c++11 -O3 -march=native -I include
LDFLAGS = -lm

N ?= 20
K ?= 3
p ?= 0.000000
s ?= 42
N_ = $(shell printf %05d $(N))
K_ = $(shell printf %04d $(K))
p_ = $(shell echo $(p) | sed 's/\./_/')
s_ = $(shell printf %d $(s))
h = "$(N_)-$(K_)-$(p_)-seed_$(s_).h"
prog = "sim-$(N_)-$(K_)-$(p_)-graphseed_$(s_)"

all: $(prog) compile_commands.json

$(h):
	./build_header.py $(N) $(K) $(p) $(s)

# make objects. '$@' = left of ':', '$^' = first item on left of ':'
$(ODIR)/%.o: %.cpp $(h)
	$(CC) -include $(h) -c -o $@ $< $(CFLAGS)

# link objects into executable '$^' = right side of ':'
$(prog): $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

compile_commands.json: Makefile
	bear $(MAKE)

.PHONY: clean
clean:
	rm -f $(prog) $(ODIR)/*.o
