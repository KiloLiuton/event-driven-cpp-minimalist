ODIR = obj
EXECNAME = simulation

_OBJ = main.o
OBJ = $(patsubst %, $(ODIR)/%, $(_OBJ))

CC = g++
CFLAGS = -Wall -std=c++11 -O3 -march=native -I include
LDFLAGS = -lm

a?=1.6
header?="00020-0003-0_0-seed_42.h"
# make objects. '$@' = left of ':', '$^' = first item on left of ':'
$(ODIR)/%.o: %.cpp
	$(CC) -D COUPLING=$(a) -include $(header) -c -o $@ $< $(CFLAGS)

# link objects into executable '$^' = right side of ':'
$(EXECNAME)$(a): $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

.PHONY: clean
clean:
	rm -f $(EXECNAME) $(ODIR)/*.o
