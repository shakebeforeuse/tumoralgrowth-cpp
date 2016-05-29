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
	  perThreadDomainBegin_(nullptr),
	  perThreadDomainEnd_(nullptr),
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
	
	domainBegin_[0] = size_;
	domainBegin_[1] = size_;
	domainEnd_[0]   = 0;
	domainEnd_[1]   = 0;
	
	
	//Alloc per-thread domain
	perThreadDomainBegin_ = new volatile int*[1];
	perThreadDomainEnd_   = new volatile int*[1];
	
	perThreadDomainBegin_[0] = new int[1];
	perThreadDomainEnd_[0]   = new int[1];
	
	//Initialize
	perThreadDomainBegin_[0][0] = size;
	perThreadDomainBegin_[0][1] = size;
	perThreadDomainEnd_[0][0]   = 0;
	perThreadDomainEnd_[0][1]   = 0;
}

void TumorAutomaton::threads(int n)
{
	//Set to sequential if n < 1
	if (n < 1)
		n = 1;
	
	if (threads_ != n)
	{
		//Remove per-thread domains
		for (int i = 0; i < threads_; ++i)
		{
			delete[] perThreadDomainBegin_[i];
			delete[] perThreadDomainEnd_[i];
		}
		
		delete[] perThreadDomainBegin_;
		delete[] perThreadDomainEnd_;
		
		//Set parameters
		threads_ = n;
		tasks_   = new std::thread[n];
		
		//Reconstruct barrier
		delete barrier_;
		barrier_ = new CyclicBarrier(n + 1);
		
		
		//Create domains for n threads. Allocate.
		perThreadDomainBegin_ = new volatile int*[n];
		perThreadDomainEnd_   = new volatile int*[n];
		
		for (int i = 0; i < n; ++i)
		{
			perThreadDomainBegin_[i] = new int[2];
			perThreadDomainEnd_[i]   = new int[2];
			
			//Initialize
			perThreadDomainBegin_[i][0] = size_;
			perThreadDomainBegin_[i][1] = size_;
			
			perThreadDomainEnd_[i][0] = 0;
			perThreadDomainEnd_[i][1] = 0;
		}
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
			
			domainBegin_[0] = std::min(domainBegin_[0], perThreadDomainBegin_[0][0]);
			domainBegin_[1] = std::min(domainBegin_[1], perThreadDomainBegin_[0][1]);
			
			domainEnd_[0] = std::max(domainEnd_[0], perThreadDomainEnd_[0][0]);
			domainEnd_[1] = std::max(domainEnd_[1], perThreadDomainEnd_[0][1]);
			
			//Change iteration direction, to avoid distortion
			if (it_ == 0)
				for (int i = domainBegin_[0]; i < domainEnd_[0]; ++i)
					for (int j = domainBegin_[1]; j < domainEnd_[1]; ++j)
						updateCell(i, j, 0);
			else
				for (int i = domainEnd_[0] - 1; i >= domainBegin_[0]; --i)
					for (int j = domainEnd_[1] - 1; j >= domainBegin_[1]; --j)
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
		{
			for (int i = 0; i < threads_; ++i)
			{
				domainBegin_[0] = std::min(domainBegin_[0], perThreadDomainBegin_[i][0]);
				domainBegin_[1] = std::min(domainBegin_[1], perThreadDomainBegin_[i][1]);
				
				domainEnd_[0] = std::max(domainEnd_[0], perThreadDomainEnd_[i][0]);
				domainEnd_[1] = std::max(domainEnd_[1], perThreadDomainEnd_[i][1]);
			}
			
			//Sync
			barrier_->await();
		}
		
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
		int delta = (domainEnd_[0] - domainBegin_[0]) / threads_;
		
		int startX = domainBegin_[0] + index       * delta;			
		int endX   = domainBegin_[0] + (index + 1) * delta;
		
		if (index + 1 == threads_)
			endX = domainEnd_[0];
		
		//Not the same than using domain{Begin|End}_ in the loop.
		//This avoid looping through just-added cells
		int startY = domainBegin_[1];
		int endY   = domainEnd_[1];
		
		//Change iteration direction, to avoid distortion
		if (it_ == 0)
			for (int i = startX; i < endX; ++i)
				for (int j = startY; j < endY; ++j)
					updateCell(i, j, index);
		else
			for (int i = endX - 1; i >= startX; --i)
				for (int j = endY - 1; j >= startY; --j)
					updateCell(i, j, index);
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
	
	//Reset domain
	domainBegin_[0] = size_;
	domainBegin_[1] = size_;
	domainEnd_[0]   = 0;
	domainEnd_[1]   = 0;
	
	//Reset per-thread domain
	for (int i = 0; i < threads_; ++i)
	{
		//Initialize
		perThreadDomainBegin_[i][0] = size_;
		perThreadDomainBegin_[i][1] = size_;
		
		perThreadDomainEnd_[i][0] = 0;
		perThreadDomainEnd_[i][1] = 0;
	}
}

//Intended to be used only externally (public)
void TumorAutomaton::cellState(int x, int y, int v)
{
	//Check if it is within bounds. Do nothing if not.
	if (0 <= x && x < size_ && 0 <= y && y < size_)
	{
		if (domainBegin_[0] > x)
			domainBegin_[0] = std::max(x, 0);
		if (domainBegin_[1] > y)
			domainBegin_[1] = std::max(y, 0);
			
		if (domainEnd_[0] <= x)
			domainEnd_[0] = std::min(x + 1, size_);
		if (domainEnd_[1] <= y)
			domainEnd_[1] = std::min(y + 1, size_);
		
		tissue_[x][y] = v;
	}
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
					
					lock_.lock();
					
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
									
									//Update per-thread domain
									if (perThreadDomainBegin_[index][0] >= x + i)
										perThreadDomainBegin_[index][0] = std::max(x + i, 0);
									if (perThreadDomainBegin_[index][1] >= y + j)
										perThreadDomainBegin_[index][1] = std::max(y + j, 0);
										
									if (perThreadDomainEnd_[index][0] <= x + i)
										perThreadDomainEnd_[index][0] = std::min(x + i + 1, size_);
									if (perThreadDomainEnd_[index][1] <= y + j)
										perThreadDomainEnd_[index][1] = std::min(y + j + 1, size_);
									
									//Stop iteration (position already chosen!)
									continueIt = false;
								}
					}
					
					lock_.unlock();
				}
			}
		}
		else
		{
			lock_.lock();
			
			//If the cell does not survive
			//Kill the cell
			tissue_[x][y] = DEAD;
			
			//Mark DORMANT neighbours as ALIVE, to be processed
			awakeNeighbourhood(x, y);
			
			lock_.unlock();
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
	
	delete barrier_;
	
	for (int i = 0; i < threads_; ++i)
	{
		delete[] perThreadDomainBegin_[i];
		delete[] perThreadDomainEnd_[i];
	}
	
	delete[] perThreadDomainBegin_;
	delete[] perThreadDomainEnd_;
}
