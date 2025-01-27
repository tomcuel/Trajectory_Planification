CC=g++ -std=c++17
CGFLAGS=-Wall -Wfatal-errors 

all : test_point.x

test_point.x : test_point.o ../../Src/point_class.o 
	$(CC) $(CGFLAGS) -o $@ $^

test_point.o: test_point.cpp
	$(CC) $(CGFLAGS) -o $@ -c $< 

../../Src/point_class.o: ../../Src/point_class.cpp
	$(CC) $(CGFLAGS) -o $@ -c $<

clean:
	rm -f test_point.o
	rm -f ../../Src/point_class.o

realclean: clean
	rm -f *.x