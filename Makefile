all: packed

packed: packed.exe

packed.exe: packed.cpp
	x86_64-w64-mingw32-g++ -o packed.exe -static -static-libgcc -static-libstdc++ packed.cpp

clean:
	rm packed.exe