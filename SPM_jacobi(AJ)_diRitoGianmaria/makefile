CXX = g++
CXXFLAGS = std=c++20 -I /usr/local/include/fastflow/ -I /usr/local/include/grppi/include/
LDFLAGS = -pthread

jacobi: jacobi.cpp
	   g++ -std=c++20 -O3 -pthread -I /usr/local/include/fastflow/ -I /usr/local/include/grppi/include/ jacobi.cpp -o jacobi