all: linux windows

linux: main.cpp src/neural.cpp src/checker.cpp src/progress.cpp includes/neural.hpp includes/parser.hpp includes/tensor.hpp includes/checker.hpp includes/progress.hpp
	g++ main.cpp src/neural.cpp src/checker.cpp src/progress.cpp -Iincludes -fopenmp -o lin

windows: main.cpp src/neural.cpp src/checker.cpp src/progress.cpp includes/neural.hpp includes/parser.hpp includes/tensor.hpp includes/checker.hpp includes/progress.hpp
	x86_64-w64-mingw32-g++ -static -static-libgcc -static-libstdc++ main.cpp src/neural.cpp src/checker.cpp src/progress.cpp -Iincludes -fopenmp -o win.exe

clean:
	rm win.exe lin