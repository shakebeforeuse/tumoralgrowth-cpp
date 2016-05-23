all: time

time: *.cpp
	g++ -O3 -std=c++11 -pthread -Wall -pedantic $^ -o time

clean:
	@rm -f *.o time

.PHONY: clean
