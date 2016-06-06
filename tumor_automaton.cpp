#include <atomic>
#include <random>
#include <thread>
#include "cyclic_barrier.h"
#include "tumor_automaton.h"

//Not a member function!
float rand_float()
{
	//Needed to seed the PRNG
	static thread_local std::random_device rd;
	
	//Keep one per thread, there is a linear cost associated.
	static thread_local std::mt19937 re(rd());
	
	//No perfomance increase if saved (it's inline)
    std::uniform_real_distribution<float> random(0, 1);
    
    return random(re);
}


thread_local char TumorAutomaton::it_ = 1;
thread_local char TumorAutomaton::prev_it_;

TumorAutomaton::TumorAutomaton(int size)
	: ps(.99),
	  pp(.8),
	  pm(.2),
	  np(5),
	  rho(2),
	  size_(size),
	  threads_(1),
	  tasks_(nullptr),
	  barrier_(nullptr)
{
	//Alloc CA
	tissue_     = new volatile int*[size_];
	ph_         = new volatile int*[size_];
	rhos_       = new volatile int*[size_];
	generation_ = new volatile char*[size_];
	
	for (int i = 0; i < size_; ++i)
	{
		tissue_[i]     = new int[size_];
		ph_[i]         = new int[size_];
		rhos_[i]       = new int[size_];
		generation_[i] = new char[size_];
		
		//Initialize
		for (int j = 0; j < size_; ++j)
		{
			tissue_[i][j]     = 0;
			ph_[i][j]         = 0;
			rhos_[i][j]       = 0;
			generation_[i][j] = 0;
		}
	}
}

void TumorAutomaton::threads(int n)
{
	//Set to sequential if n < 1
	if (n < 1)
		n = 1;
	
	if (threads_ != n)
	{
		//Set parameters
		threads_ = n;
		
		delete[] tasks_;
		tasks_   = new std::thread[n];
		
		delete locks_;
		locks_ = new std::mutex[n];
		
		//Reconstruct barrier
		delete barrier_;
		barrier_ = new CyclicBarrier(n + 1);
	}
}

void TumorAutomaton::execute(int nGenerations)
{
	//Sequential
	if (threads_ <= 1)
		//Iterate k times over the whole grid, updating cells
		for (int k = 0; k < nGenerations; ++k)
		{
			it_ = (it_ + 1) % 2;
			
			//Change iteration direction, to avoid distortion
			if (it_ == 0)
				for (int i = 0; i < size_; ++i)
					for (int j = 0; j < size_; ++j)
						updateCell(i, j, 0);
			else
				for (int i = size_ - 1; i >= 0; --i)
					for (int j = size_ - 1; j >= 0; --j)
						updateCell(i, j, 0);
		}
	else
	//Multithread
	{
		//Launch threads
		for (int i = 0; i < threads_; ++i)
			tasks_[i] = std::thread(&TumorAutomaton::operator(), this, i, nGenerations);
		
		//Compute domain limit at the start of each generation
		for (int k = 0; k < nGenerations; ++k)
			//Sync
			barrier_->await();
		
		//Join threads (in C++, either you join or detach them)
		for (int i = 0; i < threads_; ++i)
			tasks_[i].join();
	}
}

//"Runnable" method.
void TumorAutomaton::operator ()(int index, int nGenerations)
{
	for (int k = 0; k < nGenerations; ++k)
	{
		//Thread independent operation
		it_ = (it_ + prev_it_ + 1) % 2;
		
		//Await until main thread computes domain
		barrier_->await();
		
		//Compute domain for this thread
		int delta = size_ / threads_;
		
		int startX = index       * delta;			
		int endX   = (index + 1) * delta;
		
		if (index + 1 == threads_)
			endX = size_;
		
		
		//Change iteration direction, to avoid distortion
		if (it_ == 0)
			for (int i = startX; i < endX; ++i)
			{
				if (index != 0 && i < startX + 2)
					locks_[index - 1].lock();
						
				if (index != threads_ - 1 && i >= endX - 2)
					locks_[index].lock();
						
				for (int j = 0; j < size_; ++j)
					updateCell(i, j, index);
					
				if (index != 0 && i < startX + 2)
						locks_[index- 1].unlock();
						
				if (index != threads_ - 1 && i >= endX - 2)
					locks_[index].unlock();
			}
		else
			for (int i = endX - 1; i >= startX; --i)
			{
				if (index != 0 && i < startX + 2)
						locks_[index - 1].lock();
						
				if (index != threads_ - 1 && i >= endX - 2)
					locks_[index].lock();
					
				for (int j = size_ - 1; j >= 0; --j)
				{
					updateCell(i, j, index);
				}
				
				if (index != 0 && i < startX + 2)
						locks_[index - 1].unlock();
						
				if (index != threads_ - 1 && i >= endX - 2)
					locks_[index].unlock();
			}
	}
	
	//Save state for the next call to execute (if any)
	prev_it_ = it_;
}

void TumorAutomaton::reset()
{
	//Reset CA
	for (int i = 0; i < size_; ++i)
		for (int j = 0; j < size_; ++j)
		{
			tissue_[i][j]     = 0;
			ph_[i][j]         = 0;
			rhos_[i][j]       = 0;
			generation_[i][j] = 0;
		}
}

//Intended to be used only externally (public)
void TumorAutomaton::cellState(int x, int y, int v)
{
	//Check if it is within bounds. Do nothing if not.
	if (0 <= x && x < size_ && 0 <= y && y < size_)
		tissue_[x][y] = v;
}

int TumorAutomaton::cellState(int x, int y) const
{
	//Check if it is within bounds
	if (0 <= x && x < size_ && 0 <= y && y < size_)
		return tissue_[x][y];

	//If not, return ALIVE. CA won't proliferatete or migrate if it
	//thinks the cell is alive, and it does not alter CA behavior.
	return ALIVE;
}

void TumorAutomaton::awakeNeighbourhood(int x, int y)
{
	//Awakes DORMANT cells in the (x, y) neighbourhood
	for (int i = -1; i <= 1; ++i)
		for (int j = -1; j <= 1; ++j)
			if (0 <= x+i && x+i < size_ && 0 <= y+j && y+j < size_)
				if ((i != 0 || j != 0) && cellState(x+i, y+j) == DORMANT)
				{
					tissue_[x+i][y+j]     = ALIVE;
					generation_[x+i][y+j] = (it_ + 1) % 2;
				}
}

void TumorAutomaton::updateCell(int x, int y, int index)
{	
	//Check if ALIVE and whether should be processed or not (current generation?)
	if (tissue_[x][y] != DEAD && generation_[x][y] == it_)
	{
		//Mark to be computed in the next generation
		generation_[x][y] = (it_ + 1) % 2;
		
		//Check if survives
		if (rand_float() < ps)
		{
			//If DORMANT, do nothing. There is no room in the
			//neighbourhood
			if (tissue_[x][y] != DORMANT)
			{
				//Change state to alive (just for correct color in GUI)
				tissue_[x][y] = ALIVE;
				
				//Check proliferatetion
				bool proliferate = rand_float() < pp && ++ph_[x][y] >= np;
				
				//Check whether proliferates or migrates
				if (proliferate || rand_float() < pm)
				{
					//Compute next position
					float denom =  0;
					
					int n[8];
					float p[8];
					
					//Compute no. of alive neighbours
					int count = 0;
					for (int i = -1; i <= 1; ++i)
						for (int j = -1; j <= 1; ++j)
							if (i != 0 || j != 0)
							{
								n[count] = cellState(x+i, y+j) <= DEAD ? 1:0;
								denom += n[count++];
							}
					
					
					//denom == 0 means that neighbourhood is full.
					//Mark as DORMANT and stop updating this cell.
					if (denom == 0)
						tissue_[x][y] = DORMANT;
					else
					{
						//Compute probability of selecting each
						//neighbour cell.
						p[0] = n[0]/denom;
						for (int i = 1; i < 8; ++i)
							p[i] = n[i]/denom + p[i-1];
						
						
						//Select position, randomly
						float r = rand_float();
						
						int cont = 0;
						bool continueIt = true;
						for (int i = -1; i <= 1 && continueIt; ++i)
							for (int j = -1; j <= 1 && continueIt; ++j)
								if ((i != 0 || j != 0) && r < p[cont++])
								{
									//Proliferate (or migrate) to the specified cell
									if (proliferate)
									{
										//New cell				
										tissue_[x+i][y+j] = NEW;
										
										//Set proliferation signals to 0.
										ph_[x+i][y+j] = 0;
										
										//Reset no. of proliferations remaining until death.
										rhos_[x+i][y+j] = rho;
										
										//Kill the cell if it has reached the limit.
										if (--rhos_[x][y] == 0)
										{
											tissue_[x][y] = DEAD;
											awakeNeighbourhood(x, y);
										}
									}
									else
									{
										//Move to the selected position
										tissue_[x][y]     = DEAD;
										tissue_[x+i][y+j] = MIGRATED;
										awakeNeighbourhood(x, y);
										
										//Move no. of proliferation signals
										ph_[x+i][y+j] = ph_[x][y];
										ph_[x][y]     = 0;
										
										//Move no. of proliferations remaining.
										rhos_[x+i][y+j] = rhos_[x][y];
										rhos_[x][y]     = 0;
									}
									
									//Mark to be processed in the next iteration
									generation_[x+i][y+j] = (it_ + 1) % 2;
									
									//Stop iteration (position already chosen!)
									continueIt = false;
								}
					}
				}
			}
		}
		else
		{
			//If the cell does not survive
			//Kill the cell
			tissue_[x][y] = DEAD;
			
			//Mark DORMANT neighbours as ALIVE, to be processed
			awakeNeighbourhood(x, y);
		}
	}
}

TumorAutomaton::~TumorAutomaton()
{
	for (int i = 0; i < size_; ++i)
	{
		delete[] tissue_[i];
		delete[] ph_[i];
		delete[] rhos_[i];
		delete[] generation_[i];
	}
	
	delete[] tissue_;
	delete[] ph_;
	delete[] rhos_;
	delete[] generation_;
	delete[] tasks_;
	delete[] locks_;
	
	delete barrier_;
}
