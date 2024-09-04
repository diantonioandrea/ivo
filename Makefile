.PHONY: all lib distclean

CXXFLAGS = -Wall -pedantic -march=native -O2 -fPIC -I./include -O2 -fno-unsafe-math-optimizations -fno-fast-math

ifeq ($(shell uname),Darwin) # Looks for Homebrew GCC 14.
ifneq ($(shell find /opt/homebrew/bin | grep g++-14),)
CXX = g++-14
else # Apple's Clang.
CXXFLAGS += -ffp-model=precise
endif
endif

# C++ standard.
ifneq ($(shell $(CXX) --version | grep "GCC 14"),) # C++26, requires GCC 14+.
CXXFLAGS += -std=c++2c
else # C++23.
CXXFLAGS += -std=c++2b
endif

# Headers, recompilation purposes.
HEADERS = ./include/Ivo/*.hpp

# Headers, Base.
HEADERS += ./include/Ivo/Base/*.hpp

# Headers, Algebra.
HEADERS += ./include/Ivo/Algebra/*.hpp
HEADERS += ./include/Ivo/Algebra/Methods/*.hpp

# Headers, Geometry21.
HEADERS += ./include/Ivo/Geometry21/*.hpp
HEADERS += ./include/Ivo/Geometry21/Methods/*.hpp

# Headers, Mesh2121.
HEADERS += ./include/Ivo/Mesh21/*.hpp

# Library.
LIBRARY_NAME = Ivo
LIBRARY = ./lib/lib$(LIBRARY_NAME).a
INCLUDE_DESTINATION = $(HOME)/include
LIB_DESTINATION = $(HOME)/lib

# Tests.
TESTS = $(subst src/,executables/,$(subst .cpp,.out,$(shell find src -name "Test_*.cpp")))

# Source.
T_OBJECTS = $(subst src/,objects/,$(subst .cpp,.o,$(shell find src -name "Test_*")))
OBJECTS = $(subst src/,objects/,$(subst .cpp,.o,$(shell find src -name "*.cpp" -not -name "Test_*")))

# Directories.
DIRECTORIES = ./output ./objects ./executables

# All.
all: $(DIRECTORIES) $(TESTS)
	@echo "Compiled everything!"

# Library.
lib: $(LIBRARY)

$(LIBRARY): $(DIRECTORIES) $(OBJECTS)
	@mkdir ./lib
	@echo "Archiving the library to $(LIBRARY)"
	@ar rcs $(LIBRARY) $(OBJECTS)

ifeq ($(shell find . -wholename $(LIBRARY)), $(LIBRARY)) # Available only after compilation.
install: # Manual spacing for consistency between platforms.
	@echo "Installing the library"
	@echo "    - Includes under $(INCLUDE_DESTINATION)"
	@echo "    - Lib under $(LIB_DESTINATION)"
	@echo "Remember: -I$(INCLUDE_DESTINATION) -L$(LIB_DESTINATION) -l$(LIBRARY_NAME)"

	@mkdir -p $(INCLUDE_DESTINATION)
	@mkdir -p $(LIB_DESTINATION)

	@cp -r ./include/* $(INCLUDE_DESTINATION)/
	@cp $(LIBRARY) $(LIB_DESTINATION)/
endif

# Tests.
$(TESTS): executables/Test_%.out: objects/Test_%.o $(OBJECTS) 
	@if [ "$(LDFLAGS) $(LDLIBS)" = " " ]; then echo "Linking to $@"; else echo "Linking to $@ with: $(LDFLAGS) $(LDLIBS)"; fi
	@$(CXX) $(LDFLAGS) $(LDLIBS) $^ -o $@

$(T_OBJECTS): objects/%.o: src/%.cpp $(HEADERS)
	@echo "Compiling $< using $(CXX) with: $(CXXFLAGS) $(CPPFLAGS)"
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

# Objects.
$(OBJECTS): objects/%.o: src/%.cpp $(HEADERS)
	@echo "Compiling $< using $(CXX) with: $(CXXFLAGS) $(CPPFLAGS)"
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

# Directories.
$(DIRECTORIES):
	@mkdir -p $(DIRECTORIES)

# Clean.
distclean:
	@echo "Cleaning the repo."
	@$(RM) -r $(DIRECTORIES)
	@$(RM) -r ./lib
