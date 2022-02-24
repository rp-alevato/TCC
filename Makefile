# Separar em mais jobs, assim é possível usar mais núcleos da máquina e acelerar a compilação
MAKEFLAGS += -j4

# Directories
DIR_SOURCE:=src
DIR_BUILD:=build
DIR_EIGEN:=lib/eigen
DIRS_SOURCE:=$(shell find $(DIR_SOURCE) -type d)

# Files
SOURCES:=$(wildcard $(patsubst %,%/*.cpp, $(DIRS_SOURCE)))
OBJECTS:=$(patsubst $(DIR_SOURCE)/%,$(DIR_BUILD)/%, $(SOURCES:.cpp=.o))

# Dependencies
DEPS_LIBDOA:=$(filter $(DIR_BUILD)/libdoa/%, $(OBJECTS))
DEPS_MISC:=$(filter $(DIR_BUILD)/misc/%, $(OBJECTS))
DEPS_SAVE_SPECTRUM:=$(DIR_BUILD)/research/save_spectrum.o $(DEPS_LIBDOA) $(DEPS_MISC)
DEPS_SAVE_MUSIC_RESULTS_ANGLES:=$(DIR_BUILD)/research/save_music_result_angles.o $(DEPS_LIBDOA) $(DEPS_MISC)
DEPS_HP_ANALYSIS:=$(DIR_BUILD)/research/hp_analysis.o $(DEPS_LIBDOA) $(DEPS_MISC)
DEPS_CPP_TESTS:=$(DIR_BUILD)/cpp_tests/tests.o $(DEPS_LIBDOA) $(DEPS_MISC)
DEPS_TEMP:=$(DIR_BUILD)/research/temp_tests.o $(DEPS_LIBDOA) $(DEPS_MISC)

# Compiler and its flags
CC=g++-11
CXXFLAGS:=-W -Wall -Wno-format -pedantic -g -I $(DIR_EIGEN) -I $(DIR_SOURCE) -MMD -march=native -std=c++17 -Ofast -fopenmp

# Dependencies
DEPS:=$(OBJECTS:.o=.d)
-include $(DEPS)

all: save_spectrum save_music_result_angles hp_analysis temp

save_spectrum: $(DEPS_SAVE_SPECTRUM)
	$(CC) $(CXXFLAGS) $^ -o $@
	@echo "Done!\n"

save_music_result_angles: $(DEPS_SAVE_MUSIC_RESULTS_ANGLES)
	$(CC) $(CXXFLAGS) $^ -o $@
	@echo "Done!\n"

hp_analysis: $(DEPS_HP_ANALYSIS)
	$(CC) $(CXXFLAGS) $^ -o $@
	@echo "Done!\n"

cpp_tests: $(DEPS_CPP_TESTS)
	$(CC) $(CXXFLAGS) $^ -o $@
	@echo "Done!\n"

temp: $(DEPS_TEMP)
	$(CC) $(CXXFLAGS) $^ -o $@
	@echo "Done!\n"

$(DIR_BUILD)/%.o: $(DIR_SOURCE)/%.cpp
	@echo "\nCompilando $<:"
	@mkdir -p $(@D)
	$(CC) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf build/*
	-$(RM) save_spectrum save_music_result_angles hp_analysis temp cpp_tests
