# Compiler and compiler flags
CXX = g++
CXXFLAGS =  -g -Wall -O0 -std=c++17 -I./include -I./linearalgebra -I/usr/include/hdf5/serial/

# Directories
SRCDIR = src
OBJDIR = object

# Library directories
LIB_DIRS = -L/usr/lib/x86_64-linux-gnu/hdf5/serial/
# Libraries to link
LIBS = -lhdf5 -lhdf5_cpp

# Source files
SOURCES := $(wildcard $(SRCDIR)/*.cpp)

# Object files
OBJECTS := $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SOURCES))

# Executable name
EXECUTABLE = pic

# Default rule
all: welcome $(OBJDIR) compiling_field compiling_linalg $(EXECUTABLE)
	@echo "Compilation complete. Run './$(EXECUTABLE) input.ini' to execute."

# Welcome message
welcome:
	@echo "2D Electrostatic  PIC (Particle-in-Cell) Simulation Program Compilation!"

# Create the object directory
$(OBJDIR):
	@mkdir -p $(OBJDIR)

# Compile each source file into object files
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@echo "Compiling $< ..."
	@$(CXX) $(CXXFLAGS) -c $< -o $@


# Linking object files to create executable
$(EXECUTABLE): $(OBJECTS)
	@echo "Linking object files to create $(EXECUTABLE) ..."
	@$(CXX) $(CXXFLAGS) -o $@ $^ $(LIB_DIRS) $(LIBS)

# Clean rule to remove object files and executable
clean:
	@rm -f $(OBJDIR)/*.o $(EXECUTABLE)
	@echo "Cleanup complete."

# Phony targets
.PHONY: all welcome compiling_field compiling_linalg clean
