# Compiler and flags
CXX      := g++
CXXFLAGS := -std=c++17 -Wall -Wextra -O2 -Iinclude

# Directories
SRCDIR   := src
BUILDDIR := build
BINDIR   := bin

# Target
TARGET   := $(BINDIR)/Euler

# Source and object files
SOURCES  := $(wildcard $(SRCDIR)/*.cpp)
OBJECTS  := $(patsubst $(SRCDIR)/%.cpp,$(BUILDDIR)/%.o,$(SOURCES))

# Default target
all: $(TARGET)

# Link step
$(TARGET): $(OBJECTS) | $(BINDIR)
	$(CXX) $(OBJECTS) -o $@

# Compile step: src/*.cpp -> build/*.o
$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp | $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Ensure build/ and bin/ exist
$(BUILDDIR) $(BINDIR):
	mkdir -p $@


# Clean object files and binary keeping directories
clean:
	rm -f $(BUILDDIR)/*.o $(TARGET)

# cleanup
distclean: clean
	rm -rf $(BINDIR)

.PHONY: all clean distclean run
