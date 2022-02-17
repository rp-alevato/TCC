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
DEPS_SIM:=$(filter $(DIR_BUILD)/simulation/%, $(OBJECTS)) $(DEPS_LIBDOA) $(DEPS_MISC)
DEPS_LIVEAPP:=$(filter $(DIR_BUILD)/live_app/%, $(OBJECTS)) $(DEPS_LIBDOA)
DEPS_TESTS:=$(filter $(DIR_BUILD)/cpp_tests/%, $(OBJECTS)) $(DEPS_LIBDOA)

# Compiler and its flags
CC=g++
CXXFLAGS:=-W -Wall -Wno-format -pedantic -g -I $(DIR_EIGEN) -I $(DIR_SOURCE) -MMD -march=native -std=c++17 -Ofast

# Dependencies
DEPS:=$(OBJECTS:.o=.d)
-include $(DEPS)

all: sim

run: sim
	@echo "Running executable sim:\n"
	@./sim

sim: $(DEPS_SIM)
	$(CC) $^ -o $@
	@echo "Done!\n"

tests: $(DEPS_TESTS)
	$(CC) $^ -o $@
	@echo "Done!\n"

liveapp:
	$(CC) $^ -o $@
	@echo "Done!\n"

$(DIR_BUILD)/%.o: $(DIR_SOURCE)/%.cpp
	@echo "\nCompilando $<:"
	@mkdir -p $(@D)
	$(CC) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf build/*
	-$(RM) sim liveapp
