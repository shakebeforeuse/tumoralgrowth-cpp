all: time speedup

time: cyclic_barrier.cpp tumor_automaton.cpp time.cpp
	g++ -O3 -std=c++11 -pthread -Wall -pedantic $^ -o time

speedup: cyclic_barrier.cpp tumor_automaton.cpp speedup.cpp
	g++ -O3 -std=c++11 -pthread -Wall -pedantic $^ -o speedup

clean:
	@rm -f *.o time speedup

.PHONY: clean
