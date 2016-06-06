#ifndef TUMORAUTOMATON_H_
#define TUMORAUTOMATON_H_

#include <atomic>
#include <random>
#include <thread>
#include "cyclic_barrier.h"

class TumorAutomaton
{
	public:
		//Cell states
		static const int DEAD     = 0;
		static const int DORMANT  = 1;
		static const int ALIVE    = 2;
		static const int NEW      = 3;
		static const int MIGRATED = 4;
		
		//Tumoral growth simulation parameters
		double ps;
		double pp;
		double pm;
		int    np;
		int    rho;
		
		//API
		explicit TumorAutomaton(int);
		
		void threads(int);
		void execute(int);
		void operator()(int, int);
		void reset();
		
		int cellState(int, int) const;
		void cellState(int, int, int);
		
		~TumorAutomaton();
		

	private:
		//CA
		volatile int**  tissue_;
		volatile int**  rhos_;
		volatile int**  ph_;
		volatile char** generation_;
		
		static thread_local char it_;
		static thread_local char prev_it_;
		int    size_;
		
		//Paralelism
		//Number of threads we will have and array of tasks.
		int threads_;
		std::thread* tasks_;
		
		//Synchronization
		CyclicBarrier* barrier_;
		std::mutex* locks_;
		
		
		void awakeNeighbourhood(int, int);
		void updateCell(int, int);
		void updateCell(int, int, int);
};

#endif
