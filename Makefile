.PHONY: all distclean
CXXFLAGS = -Wall -pedantic -Wno-pre-c++2b-compat -std=c++2b -march=native -O2 -fPIC -I./include -O2 -fno-unsafe-math-optimizations -fno-fast-math

ifeq ($(shell uname),Darwin) # Apple's clang.
CXXFLAGS += -ffp-model=precise
endif

# Headers, recompilation purposes.
HEADERS = ./include/*.hpp

# Base.
HEADERS += ./include/Base/*.hpp

# Algebra.
HEADERS += ./include/Algebra/*.hpp
HEADERS += ./include/Algebra/Methods/*.hpp

# Tests.
EXECS = $(subst test/,executables/,$(subst .cpp,.out,$(shell find test -name "*.cpp")))
OBJECTS = $(subst test/,objects/,$(subst .cpp,.o,$(shell find test -name "*.cpp")))

# Directories.
DIRECTORIES = ./output ./objects ./executables

# All.
all: $(DIRECTORIES) $(EXECS)
	@echo "Compiled everything!"

$(EXECS): executables/%.out: objects/%.o $(OBJECTS) 
	@if [ "$(LDFLAGS) $(LDLIBS)" = " " ]; then echo "Linking $(subst objects/,,$<) and base objects to $@"; else echo "Linking $(subst objects/,,$<) and base objects to $@ with: $(LDFLAGS) $(LDLIBS)"; fi
	@$(CXX) $(LDFLAGS) $(LDLIBS) $^ -o $@

$(OBJECTS): objects/%.o: test/%.cpp $(HEADERS)
	@echo "Compiling $< using $(CXX) with: $(CXXFLAGS) $(CPPFLAGS)"
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

# Directories.
$(DIRECTORIES):
	@mkdir -p $(DIRECTORIES)

# Clean.
distclean:
	@echo "Cleaning the repo."
	@$(RM) -r $(DIRECTORIES)
