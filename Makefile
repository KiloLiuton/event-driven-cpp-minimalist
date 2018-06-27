ODIR = obj

_OBJ = main.o
OBJ = $(patsubst %, $(ODIR)/%, $(_OBJ))

CC = g++
CFLAGS = -Wall -std=c++11 -O3 -march=native -I include
LDFLAGS = -lm

a ?= 1.6
N ?= 20
K ?= 3
p ?= 0.0
s ?= 42
a_ = $(shell echo $(a) | sed 's/\./_/')
N_ = $(shell printf %05d $(N))
K_ = $(shell printf %04d $(K))
p_ = $(shell echo $(p) | sed 's/\./_/')
s_ = $(shell printf %d $(s))
h = "$(N_)-$(K_)-$(p_)-seed_$(s_).h"
prog = "sim-$(N_)-$(K_)-$(p_)-$(a_)-graphseed_$(s_)"

all: $(prog)

$(h):
	./build_header.py $(N) $(K) $(p) $(s)

# make objects. '$@' = left of ':', '$^' = first item on left of ':'
$(ODIR)/%.o: %.cpp $(h)
	$(CC) -D COUPLING=$(a) -include $(h) -c -o $@ $< $(CFLAGS)

# link objects into executable '$^' = right side of ':'
$(prog): $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

.PHONY: clean
clean:
	rm -f $(prog) $(ODIR)/*.o
