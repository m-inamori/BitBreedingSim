CXX = g++
CXXFLAGS = -std=gnu++17 -I"/usr/share/R/include" -I"/home/inamori/R/x86_64-pc-linux-gnu-library/4.4/Rcpp/include" -O3 -fpic
LDFLAGS = -shared -L/usr/lib/R/lib -lR

BUILD_DIR = build/release
SRC_DIR = src
TARGET = $(BUILD_DIR)/BitBreedingSim.so

SOURCES = $(wildcard src/*.cpp)
OBJECTS = $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SOURCES))

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) $(LDFLAGS) -o $(TARGET)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
    mkdir -p $(BUILD_DIR)
    $(CXX) $(CXXFLAGS) -c $< -o $@

install:
	mkdir -p libs
	cp $(BUILD_DIR)/BitBreedingSim.so libs/
	R CMD INSTALL .

clean:
    rm -rf $(BUILD_DIR)
