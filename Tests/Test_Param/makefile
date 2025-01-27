CC=g++ -std=c++17
CGFLAGS=-Wall -Wfatal-errors 

all : test_param.x

test_param.x : test_param.o
	$(CC) $(CGFLAGS) -o $@ $^

test_param.o: test_param.cpp
	$(CC) $(CGFLAGS) -o $@ -c $< 

clean:
	rm -f test_param.o

realclean: clean
	rm -f *.x