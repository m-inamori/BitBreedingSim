CXX = g++
include build/Makevars

BUILD_DIR = build/release
SRC_DIR = src
PKG_LIBS = -L$(LDFLAGS) -lBitBreedingSim
TARGET = $(BUILD_DIR)/BitBreedingSim.so

SOURCES = $(wildcard src/*.cpp)
OBJECTS = $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SOURCES))

CXXFLAGS = -std=gnu++17 $(PKG_CPPFLAGS) -O3 -fpic

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) $(LDFLAGS) -o $(TARGET)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

install:
	R CMD INSTALL .

clean:
	rm -rf $(BUILD_DIR)
