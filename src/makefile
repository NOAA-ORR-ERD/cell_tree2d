all: cell_tree2d

cell_tree2d: testmain.o cell_tree2d.o
	g++ -g *.o -o cell_tree2dtest.exe

cell_tree2d.o: cell_tree2d.cpp
	g++ -g -c cell_tree2d.cpp

testmain.o: testmain.cpp
	g++ -g -std=c++11  -c testmain.cpp

clean:
	rm *o cell_tree2dtest.exe
	
