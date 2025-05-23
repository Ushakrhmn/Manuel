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

BIN =  probability_NOvA_NH_dcp-90_newphysics probability_NOvA_NH_dcp-90_SI probability_NOvA_NH_dcp90_newphysics probability_NOvA_NH_dcp90_SI probability_T2K_NH_dcp-90_newphysics probability_T2K_NH_dcp-90_SI probability_T2K_NH_dcp90_newphysics probability_T2K_NH_dcp90_SI
OBJ = probability_NOvA_NH_dcp-90_newphysics.o probability_NOvA_NH_dcp-90_newphysics.o probability_NOvA_NH_dcp90_newphysics.o probability_NOvA_NH_dcp90_SI.o probability_T2K_NH_dcp-90_newphysics.o probability_T2K_NH_dcp-90_SI.o probability_T2K_NH_dcp90_newphysics.o probability_T2K_NH_dcp90_SI.o


probability_NOvA_NH_dcp-90_newphysics: probability_NOvA_NH_dcp-90_newphysics.o myio.o
	gcc probability_NOvA_NH_dcp-90_newphysics.o myio.o -o probability_NOvA_NH_dcp-90_newphysics $(LDFLAGS) $(local_LDFLAGS) $(PYTHON_LIB)

probability_NOvA_NH_dcp90_newphysics: probability_NOvA_NH_dcp90_newphysics.o myio.o
	gcc probability_NOvA_NH_dcp90_newphysics.o myio.o -o probability_NOvA_NH_dcp90_newphysics $(LDFLAGS) $(local_LDFLAGS) $(PYTHON_LIB)
probability_NOvA_NH_dcp-90_SI: probability_NOvA_NH_dcp-90_SI.o myio.o
	gcc probability_NOvA_NH_dcp-90_SI.o myio.o -o probability_NOvA_NH_dcp-90_SI $(LDFLAGS) $(local_LDFLAGS) $(PYTHON_LIB)

probability_NOvA_NH_dcp90_SI: probability_NOvA_NH_dcp90_SI.o myio.o
	gcc probability_NOvA_NH_dcp90_SI.o myio.o -o probability_NOvA_NH_dcp90_SI $(LDFLAGS) $(local_LDFLAGS) $(PYTHON_LIB)

probability_T2K_NH_dcp-90_newphysics: probability_T2K_NH_dcp-90_newphysics.o myio.o
	gcc probability_T2K_NH_dcp-90_newphysics.o myio.o -o probability_T2K_NH_dcp-90_newphysics $(LDFLAGS) $(local_LDFLAGS) $(PYTHON_LIB)

probability_T2K_NH_dcp-90_SI: probability_T2K_NH_dcp-90_SI.o myio.o
	gcc probability_T2K_NH_dcp-90_SI.o myio.o -o probability_T2K_NH_dcp-90_SI $(LDFLAGS) $(local_LDFLAGS) $(PYTHON_LIB)

probability_T2K_NH_dcp90_newphysics: probability_T2K_NH_dcp90_newphysics.o myio.o
	gcc probability_T2K_NH_dcp90_newphysics.o myio.o -o probability_T2K_NH_dcp90_newphysics $(LDFLAGS) $(local_LDFLAGS) $(PYTHON_LIB)

probability_T2K_NH_dcp90_SI: probability_T2K_NH_dcp90_SI.o myio.o
	gcc probability_T2K_NH_dcp90_SI.o myio.o -o probability_T2K_NH_dcp90_SI $(LDFLAGS) $(local_LDFLAGS) $(PYTHON_LIB)

%.o : %.c
	gcc -std=c99 $(CFLAGS) $(local_CFLAGS) -c $< $(INCFLAGS) $(PYTHON_INC)

.PHONY: clean
clean:
	rm -f $(BIN) $(OBJ)
