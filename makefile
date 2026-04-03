CXX = g++
CXXFLAGS = -O3 -march=native -ffast-math -std=c++17

TARGET = MultiStrainSIRonNet_2patches.exe
SRC = MultiStrainSIRonNet.cpp

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET)

clean:
	rm -f $(TARGET)