.PHONY: all lib distclean

ifeq ($(shell uname),Darwin) # Forces GCC (g++-14, homebrew) instead of Clang.
CXX = g++-14
endif

CXXFLAGS += -Wall -std=c++23 -pedantic -I./include -fPIC -march=native -Ofast -fopenmp
LDLIBS += -lgomp

# Headers, recompilation purposes.
HEADERS = ./include/*.hpp
HEADERS = ./include/Ivo/*.hpp

# Headers, Base.
HEADERS += ./include/Ivo/Base/*.hpp

# Headers, Algebra.
HEADERS += ./include/Ivo/Algebra/*.hpp
HEADERS += ./include/Ivo/Algebra/Methods/*.hpp

# Headers, Geometry21.
HEADERS += ./include/Ivo/Geometry21/*.hpp
HEADERS += ./include/Ivo/Geometry21/Methods/*.hpp

# Headers, Mesh21.
HEADERS += ./include/Ivo/Mesh21/*.hpp

# Headers, Fem.
HEADERS += ./include/Ivo/Fem/*.hpp

# Headers, Problem.
HEADERS += ./include/Ivo/Problem/*.hpp

# Headers, Tests.
HEADERS += ./src/include/*.hpp

# Library.
LIBRARY_NAME = Ivo
LIBRARY = ./lib/lib$(LIBRARY_NAME).a
INCLUDE_DESTINATION = $(HOME)/include
LIB_DESTINATION = $(HOME)/lib

# Tests and scripts.
TESTS = $(subst src/,executables/,$(subst .cpp,.out,$(shell find src -name "Test_*.cpp")))
SCRIPTS = $(subst src/,executables/,$(subst .cpp,.out,$(shell find src -name "Script_*.cpp")))

# Source.
T_OBJECTS = $(subst src/,objects/,$(subst .cpp,.o,$(shell find src -name "Test_*")))
S_OBJECTS = $(subst src/,objects/,$(subst .cpp,.o,$(shell find src -name "Script_*")))

OBJECTS = $(subst src/,objects/,$(subst .cpp,.o,$(shell find src -name "Ivo_*")))

# Directories.
DIRECTORIES = ./output ./objects ./executables

# All.
all: $(DIRECTORIES) $(TESTS) $(SCRIPTS)
	@echo "Compiled everything!"

# Tests.
tests: $(DIRECTORIES) $(TESTS)
	@echo "Compiled tests!"

# Scripts.
scripts: $(DIRECTORIES) $(SCRIPTS)
	@echo "Compiled scripts!"

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

# Tests and scripts.
$(TESTS): executables/Test_%.out: objects/Test_%.o $(OBJECTS) 
	@if [ "$(LDFLAGS) $(LDLIBS)" = " " ]; then echo "Linking to $@"; else echo "Linking to $@ with: $(LDFLAGS) $(LDLIBS)"; fi
	@$(CXX) $(LDFLAGS) $(LDLIBS) $^ -o $@

$(SCRIPTS): executables/Script_%.out: objects/Script_%.o $(OBJECTS) 
	@if [ "$(LDFLAGS) $(LDLIBS)" = " " ]; then echo "Linking to $@"; else echo "Linking to $@ with: $(LDFLAGS) $(LDLIBS)"; fi
	@$(CXX) $(LDFLAGS) $(LDLIBS) $^ -o $@

$(T_OBJECTS): objects/%.o: src/%.cpp $(HEADERS)
	@echo "Compiling $< using $(CXX) with: $(CXXFLAGS) $(CPPFLAGS)"
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

$(S_OBJECTS): objects/%.o: src/%.cpp $(HEADERS)
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
