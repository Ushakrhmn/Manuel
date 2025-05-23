# Compiler and flags
CC = gcc
CFLAGS = -std=c99 -O3 -fPIC -Wall \
    -I$(shell python3 -c "from sysconfig import get_paths as gp; print(gp()['include'])") \
    -I/usr/local/include

# Link to Python3 library
PYTHON_LIB_PATH = $(shell python3 -c "from sysconfig import get_config_var as gcv; print(gcv('LIBDIR'))")
PYTHON_VERSION = $(shell python3 -c "import sys; print(sys.version_info[0])").$(shell python3 -c "import sys; print(sys.version_info[1])")
PYTHON_LIB = $(PYTHON_LIB_PATH)/libpython$(PYTHON_VERSION).so

# Explicit path to libglobes.so
LIBGLOBES_PATH = /eos/user/u/urahaman/SWAN_projects/long-baseline/globes/globes-3.2.18/lib

# If you're using NumPy 2.x, use the correct NumPy include path
NUMPY_PATH = $(shell python3 -c "import numpy; print(numpy.get_include())")

# Linker flags, include NumPy and Python, and explicitly link to libglobes
LDFLAGS = -L/usr/local/lib -L$(LIBGLOBES_PATH) \
    -L$(shell python3 -c "from sysconfig import get_paths as gp; print(gp()['stdlib'])") \
    -lglobes -lgsl -lgslcblas -lm \
    -ldl $(PYTHON_LIB)  # Include the Python library explicitly

# Optional: You may want to ensure you're linking against the correct NumPy version
NUMPY_LIB = $(shell python3 -c "import numpy; print(numpy.__file__)")
LDFLAGS += -L$(NUMPY_LIB)/../..

# Source and target settings
SRC = chi2_NH_NOvA_2024_newphysics.c
OBJ = $(SRC:.c=.o)
TARGET = chi2_NH_NOvA_2024_newphysics  # Change the target to an executable


# Default rule to build the executable
all: $(TARGET)

# Rule to create the executable from the object file
$(TARGET): $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

# Rule to compile the object file from the source
$(OBJ): $(SRC)
	$(CC) $(CFLAGS) -o $@ -c $^


