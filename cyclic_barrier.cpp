#include <condition_variable>
#include <mutex>
#include "cyclic_barrier.h"

CyclicBarrier::CyclicBarrier(unsigned initial_count)
	: backwards_(true),
	  count_(initial_count),
	  initial_count_(initial_count)
{}

void CyclicBarrier::await()
{
	std::unique_lock<std::mutex> lock_(mutex_);
	
	if (backwards_)	
		if (--count_ == 0)
		{
			all_threads_ready_.notify_all();
			backwards_ = false;
		}
		else
			//Lambda predicate to avoid spurious wake ups
			all_threads_ready_.wait(lock_, [this] {return !backwards_;});
	else
		if (++count_ == initial_count_)
		{
			all_threads_ready_.notify_all();
			backwards_ = true;
		}
		else
			//Lambda predicate to avoid spurious wake ups
			all_threads_ready_.wait(lock_, [this] {return backwards_;});
}
