# Separar em mais jobs, assim é possível usar mais núcleos da máquina e acelerar a compilação
MAKEFLAGS += -j4

# Directories
DIR_SOURCE:=src
DIR_BUILD:=build
DIR_RESEARCH:=$(DIR_BUILD)/research
DIR_EIGEN:=lib/eigen
DIRS_SOURCE:=$(shell find $(DIR_SOURCE) -type d)

# Files
SOURCES:=$(wildcard $(patsubst %,%/*.cpp, $(DIRS_SOURCE)))
OBJECTS:=$(patsubst $(DIR_SOURCE)/%,$(DIR_BUILD)/%, $(SOURCES:.cpp=.o))

# Dependencies
DEPS_AOA:=$(filter $(DIR_BUILD)/doa/%, $(OBJECTS))
DEPS_MISC:=$(filter $(DIR_BUILD)/misc/%, $(OBJECTS))
DEPS_SAVE_SPECTRUM:=$(DIR_RESEARCH)/save_spectrum.o $(DEPS_AOA) $(DEPS_MISC)
DEPS_SAVE_MUSIC_RESULTS_ANGLES:=$(DIR_RESEARCH)/save_music_result_angles.o $(DEPS_AOA) $(DEPS_MISC)
DEPS_HP_ANALYSIS:=$(DIR_RESEARCH)/hp_analysis.o $(DEPS_AOA) $(DEPS_MISC)
DEPS_PRECISION_ANALYSIS:=$(DIR_RESEARCH)/precision_analysis.o $(DEPS_AOA) $(DEPS_MISC)
DEPS_CPP_TESTS:=$(DIR_RESEARCH)/cpp_tests.o $(DEPS_AOA) $(DEPS_MISC)
DEPS_SMALL_TESTS:=$(DIR_RESEARCH)/small_tests.o $(DEPS_AOA) $(DEPS_MISC)

# Compiler and its flags
CC=g++-10
CXXFLAGS:=-W -Wall -Wno-format -pedantic -g -I $(DIR_EIGEN) -I $(DIR_SOURCE) -MMD -march=native -std=c++17 -fopenmp -O2

# Dependencies
DEPS:=$(OBJECTS:.o=.d)
-include $(DEPS)

all: save_spectrum save_music_result_angles hp_analysis precision_analysis

save_spectrum: $(DEPS_SAVE_SPECTRUM)
	$(CC) $(CXXFLAGS) $^ -o $@.exe
	@echo "Done!\n"

save_music_result_angles: $(DEPS_SAVE_MUSIC_RESULTS_ANGLES)
	$(CC) $(CXXFLAGS) $^ -o $@.exe
	@echo "Done!\n"

hp_analysis: $(DEPS_HP_ANALYSIS)
	$(CC) $(CXXFLAGS) $^ -o $@.exe
	@echo "Done!\n"

precision_analysis: $(DEPS_PRECISION_ANALYSIS)
	$(CC) $(CXXFLAGS) $^ -o $@.exe
	@echo "Done!\n"

cpp_tests: $(DEPS_CPP_TESTS)
	$(CC) $(CXXFLAGS) $^ -o $@.exe
	@echo "Done!\n"

small_tests: $(DEPS_SMALL_TESTS)
	$(CC) $(CXXFLAGS) $^ -o $@.exe
	@echo "Done!\n"

$(DIR_BUILD)/%.o: $(DIR_SOURCE)/%.cpp
	@echo "\nCompilando $<:"
	@mkdir -p $(@D)
	$(CC) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf build/*
	-$(RM) *.exe
