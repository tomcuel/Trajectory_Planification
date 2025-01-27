CC=g++ -std=c++17
CGFLAGS=-Wall -Wfatal-errors 

all : test_graph.x

test_graph.x : test_graph.o ../../Src/point_class.o ../../Src/arc_class.o ../../Src/obstacle_class.o ../../Src/graph_class.o
	$(CC) $(CGFLAGS) -o $@ $^

test_graph.o: test_graph.cpp
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
	rm -f test_graph.o
	rm -f ../../Src/point_class.o
	rm -f ../../Src/arc_class.o
	rm -f ../../Src/obstacle_class.o
	rm -f ../../Src/graph_class.o

realclean: clean
	rm -f *.x