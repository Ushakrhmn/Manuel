prefix = /usr/local
exec_prefix = /usr/local
libdir = ${exec_prefix}/lib
globesconf= $(exec_prefix)/bin/globes-config

local_CFLAGS = -g -O3

INCFLAGS:=$(shell $(globesconf) --include)
local_LDFLAGS:=$(shell $(globesconf) --libs) 
preopen_modules:=$(shell $(globesconf) --dlpreopen)
local_LTLDFLAGS:=$(shell $(globesconf) --ltlibs)
ifdef preopen_modules
predefs = -DGLB_STATIC
endif

# Python integration using venv for runtime, Homebrew for linking
PYTHON_VENV = /Users/ushak/venvs/globes
PYTHON_INC = -I/opt/homebrew/opt/python@3.13/Frameworks/Python.framework/Versions/3.13/include/python3.13
PYTHON_LIB = -L/opt/homebrew/opt/python@3.13/Frameworks/Python.framework/Versions/3.13/lib -lpython3.13 -ldl -framework CoreFoundation

BIN =  probability_NOvA_NH_dcp-90_newphysics
OBJ = probability_NOvA_NH_dcp-90_newphysics.o 

$(BIN): $(OBJ) myio.o
	gcc $(OBJ) myio.o -o $(BIN) $(LDFLAGS) $(local_LDFLAGS) $(PYTHON_LIB)

%.o : %.c
	gcc -std=c99 $(CFLAGS) $(local_CFLAGS) -c $< $(INCFLAGS) $(PYTHON_INC)

.PHONY: clean
clean:
	rm -f $(BIN) $(OBJ)