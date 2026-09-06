CXX = g++
CXXFLAGS = -O3 -ffast-math -std=c++17

TARGET = MultiStrainSIRonNet.exe
SRC = MultiStrainSIRonNet.cpp

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET)

clean:
	rm -f $(TARGET)
