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
DEPS_SIM:=$(filter $(DIR_BUILD)/sim/%, $(OBJECTS)) $(DEPS_LIBDOA)
DEPS_LIVEAPP:=$(filter $(DIR_BUILD)/live_app/%, $(OBJECTS)) $(DEPS_LIBDOA)

# Compiler and its flags
CC=g++
CXXFLAGS:=-W -Wall -Wno-format -pedantic -g -Ofast -I$(DIR_EIGEN) -I$(DIR_SOURCE) -MMD -march=native

# Dependencies
DEPS:=$(DEPS_LIBDOA:.o=.d)
include $(DEPS)

all: sim liveapp

sim: $(DEPS_SIM)
	$(CC) $^ -o $@
	@echo "Done!"

liveapp:
	$(CC) $^ -o $@
	@echo "Done!"

$(DIR_BUILD)/%.o: $(DIR_SOURCE)/%.cpp
	@echo "\nCompilando $<:"
	@mkdir -p $(@D)
	$(CC) $(CXXFLAGS) -c $< -o $@

clean:
	-$(RM) $(OBJECTS)
	-$(RM) $(DEPENDS)
	-$(RM) sim liveapp
