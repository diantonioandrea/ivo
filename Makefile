.PHONY: all distclean
CXXFLAGS = -Wall -pedantic -Wno-pre-c++2b-compat -std=c++2b -march=native -O2 -fPIC -I./include -O2 -fno-unsafe-math-optimizations -fno-fast-math

ifeq ($(shell uname),Darwin) # Apple's clang.
CXXFLAGS += -ffp-model=precise
endif

# Headers, recompilation purposes.
HEADERS = ./include/*.hpp

# Headers, Base.
HEADERS += ./include/Base/*.hpp

# Headers, Algebra.
HEADERS += ./include/Algebra/*.hpp
HEADERS += ./include/Algebra/Methods/*.hpp

# Headers, Geometry21.
HEADERS += ./include/Geometry21/*.hpp
HEADERS += ./include/Geometry21/Methods/*.hpp

# Headers, Mesh2121.
HEADERS += ./include/Mesh21/*.hpp

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

$(TESTS): executables/Test_%.out: objects/Test_%.o $(OBJECTS) 
	@if [ "$(LDFLAGS) $(LDLIBS)" = " " ]; then echo "Linking to $@"; else echo "Linking to $@ with: $(LDFLAGS) $(LDLIBS)"; fi
	@$(CXX) $(LDFLAGS) $(LDLIBS) $^ -o $@

$(T_OBJECTS): objects/%.o: src/%.cpp $(HEADERS)
	@echo "Compiling $< using $(CXX) with: $(CXXFLAGS) $(CPPFLAGS)"
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

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
