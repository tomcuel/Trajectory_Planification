CC=g++ -std=c++17
CGFLAGS=-Wall -Wfatal-errors 

all : test_dijkstra.x

test_dijkstra.x : test_dijkstra.o ../../Src/point_class.o ../../Src/arc_class.o ../../Src/obstacle_class.o ../../Src/graph_class.o
	$(CC) $(CGFLAGS) -o $@ $^

test_dijkstra.o: test_dijkstra.cpp
	$(CC) $(CGFLAGS) -o $@ -c $< 

../../Src/point_class.o: ../../Src/point_class.cpp
	$(CC) $(CGFLAGS) -o $@ -c $<

../../Src/arc_class.o: ../../Src/arc_class.cpp
	$(CC) $(CGFLAGS) -o $@ -c $<

../../Src/obstacle_class.o: ../../Src/obstacle_class.cpp
	$(CC) $(CGFLAGS) -o $@ -c $<

../../Src/graph_class.o: ../../Src/graph_class.cpp
	$(CC) $(CGFLAGS) -o $@ -c $<

clean:
	rm -f test_dijkstra.o
	rm -f ../../Src/point_class.o
	rm -f ../../Src/arc_class.o
	rm -f ../../Src/obstacle_class.o
	rm -f ../../Src/graph_class.o

realclean: clean
	rm -f *.x