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
		
		const int cellState(int, int);
		void cellState(int, int, int);
		
		~TumorAutomaton();
		

	private:
		//CA
		int**  tissue_;
		int**  rhos_;
		int**  ph_;
		char** generation_;
		char   it_;
		int    size_;
		
		//Dynamic domain
		int domainBegin_[2];
		int domainEnd_[2];
		
		//Paralelism
		//Number of threads we will have and array of tasks.
		int threads_;
		std::thread* tasks_;
		
		//Synchronization
		CyclicBarrier* barrier_;
		
		//Non-static random number generaton (to avoid thread-safety)
		std::uniform_real_distribution<float> random_;
		std::default_random_engine re_;
		
		void awakeNeighbourhood(int, int);
		void updateCell(int, int);
};

#endif
