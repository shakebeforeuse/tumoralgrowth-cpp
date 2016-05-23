#ifndef CYCLIC_BARIER_H_
#define CYCLIC_BARIER_H_

#include <condition_variable>
#include <mutex>

class CyclicBarrier
{
	public:
		explicit CyclicBarrier(unsigned);
		void await();
		
	private:
		bool backwards_;
		unsigned count_;
		const unsigned initial_count_;
		std::mutex mutex_;
		std::condition_variable all_threads_ready_;
};

#endif
