CC=g++ -std=c++17
CGFLAGS=-Wall -Wfatal-errors 

all : test_arc.x

test_arc.x : test_arc.o ../../Src/point_class.o ../../Src/arc_class.o 
	$(CC) $(CGFLAGS) -o $@ $^

test_arc.o: test_arc.cpp
	$(CC) $(CGFLAGS) -o $@ -c $< 

../../Src/point_class.o: ../../Src/point_class.cpp
	$(CC) $(CGFLAGS) -o $@ -c $<

../../Src/arc_class.o: ../../Src/arc_class.cpp
	$(CC) $(CGFLAGS) -o $@ -c $<

clean:
	rm -f test_arc.o
	rm -f ../../Src/point_class.o
	rm -f ../../Src/arc_class.o

realclean: clean
	rm -f *.x