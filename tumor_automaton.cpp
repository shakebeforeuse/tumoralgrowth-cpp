#include <atomic>
#include <random>
#include <thread>
#include "cyclic_barrier.h"
#include "tumor_automaton.h"

TumorAutomaton::TumorAutomaton(int size)
	: ps(.99),
	  pp(.8),
	  pm(.2),
	  np(5),
	  rho(2),
	  it_(0),
	  size_(size),
	  threads_(0),
	  tasks_(nullptr),
	  barrier_(nullptr)
{
	tissue_     = new int*[size_];
	ph_         = new int*[size_];
	rhos_       = new int*[size_];
	generation_ = new char*[size_];
	
	for (int i = 0; i < size_; ++i)
	{
		tissue_[i]     = new int[size_];
		ph_[i]         = new int[size_];
		rhos_[i]       = new int[size_];
		generation_[i] = new char[size_];
		
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
}

void TumorAutomaton::threads(int n)
{
	threads_ = n;
	tasks_   = new std::thread[n];
	
	delete barrier_;
	barrier_ = new CyclicBarrier(n + 1);
}

void TumorAutomaton::execute(int nGenerations)
{
	if (threads_ <= 1)
		//Iterate k times over the whole grid, updating cells
		for (int k = 0; k < nGenerations; ++k)
		{
			//Change iteration direction, to avoid distortion
			if (it_ == 0)
				for (int i = domainBegin_[0]; i < domainEnd_[0]; ++i)
					for (int j = domainBegin_[1]; j < domainEnd_[1]; ++j)
						updateCell(i, j);
			else
				for (int i = domainEnd_[0] - 1; i >= domainBegin_[0]; --i)
					for (int j = domainEnd_[1] - 1; j >= domainBegin_[1]; --j)
						updateCell(i, j);
			
			it_ = (it_ + 1) % 2;
		}
	else
	{
		for (int i = 0; i < threads_; ++i)
			tasks_[i] = std::thread(&TumorAutomaton::operator(), this, i, nGenerations);
		
		for (int k = 0; k < nGenerations; ++k)
		{
			barrier_->await();
			it_ = (it_ + 1) % 2;
		}
		
		for (int i = 0; i < threads_; ++i)
			tasks_[i].join();
	}
}

void TumorAutomaton::operator ()(int index, int nGenerations)
{
	for (int k = 0; k < nGenerations; ++k)
	{
		int delta = (domainEnd_[0] - domainBegin_[0]) / threads_;
		
		int startX = domainBegin_[0] + index       * delta;			
		int endX   = domainBegin_[0] + (index + 1) * delta;
		
		if (index + 1 == threads_)
			endX = domainEnd_[0];
		
		//Not the same than using domain{Begin|End}_ in the loop.
		//This avoid looping through just added cells
		int startY = domainBegin_[1];
		int endY   = domainEnd_[1];
		
		//Change iteration direction, to avoid distortion
		if (it_ == 0)
			for (int i = startX; i < endX; ++i)
				for (int j = startY; j < endY; ++j)
						updateCell(i, j);
		else
			for (int i = endX - 1; i >= startX; --i)
				for (int j = endY - 1; j >= startY; --j)
					updateCell(i, j);
					
		barrier_->await();
	}
}

void TumorAutomaton::reset()
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
	
	tissue_     = new int*[size_];
	ph_         = new int*[size_];
	rhos_       = new int*[size_];
	generation_ = new char*[size_];
	
	for (int i = 0; i < size_; ++i)
	{
		tissue_[i]     = new int[size_];
		ph_[i]         = new int[size_];
		rhos_[i]       = new int[size_];
		generation_[i] = new char[size_];
		
		for (int j = 0; j < size_; ++j)
		{
			tissue_[i][j]  = 0;
			ph_[i][j]      = 0;
			rhos_[i][j]    = 0;
			generation_[i] = 0;
		}
	}
	
	domainBegin_[0] = size_;
	domainBegin_[1] = size_;
	domainEnd_[0]   = 0;
	domainEnd_[1]   = 0;
}

void TumorAutomaton::cellState(int x, int y, int v)
{
	//Check if it is within bounds. Do nothing if not.
	if (0 <= x && x < size_ && 0 <= y && y < size_)
	{
		if (domainBegin_[0] > x)
			domainBegin_[0] = std::max(x - 1, 0);
		if (domainBegin_[1] > y)
			domainBegin_[1] = std::max(y - 1, 0);
			
		if (domainEnd_[0] <= x)
			domainEnd_[0] = std::min(x + 1, size_);
		if (domainEnd_[1] <= y)
			domainEnd_[1] = std::min(y + 1, size_);
		
		tissue_[x][y] = v;
	}
}

const int TumorAutomaton::cellState(int x, int y)
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
	for (int i = -1; i <= 1; ++i)
		for (int j = -1; j <= 1; ++j)
			if (0 <= x+i && x+i < size_ && 0 <= y+j && y+j < size_)
				if ((i != 0 || j != 0) && cellState(x+i, y+j) == DORMANT)
				{
					cellState(x+i, y+j, ALIVE);
					generation_[x+i][y+j] = (it_ + 1) % 2;
				}
}

void TumorAutomaton::updateCell(int x, int y)
{
	//Check if ALIVE and whether should be processed or not (current generation?)
	if (cellState(x, y) != DEAD && generation_[x][y] == it_)
	{
		//Mark to be computed in the next generation
		generation_[x][y] = (it_ + 1) % 2;
		
		//Check if survives
		if (random_(re_) < ps)
		{
			//If DORMANT, do nothing. There is no room in the
			//neighbourhood
			if (cellState(x, y) != DORMANT)
			{
				//Change state to alive (just for correct color in GUI)
				cellState(x, y, ALIVE);
				
				//Check proliferatetion
				bool proliferate = random_(re_) < pp && ++ph_[x][y] >= np;
				
				//Check whether proliferates or migrates
				if (proliferate || random_(re_) < pm)
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
								n[count] = cellState(x + i, y + j) <= DEAD ? 1:0;
								denom += n[count++];
							}
					
					
					//denom == 0 means that neighbourhood is full.
					//Mark as DORMANT and stop updating this cell.
					if (denom == 0)
						cellState(x, y, DORMANT);
					else
					{
						//Compute probability of selecting each
						//neighbour cell.
						p[0] = n[0]/denom;
						for (int i = 1; i < 8; ++i)
							p[i] = n[i]/denom + p[i-1];
						
						
						//Select position, randomly
						float r = random_(re_);
						
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
										cellState(x + i, y + j, NEW);
										
										//Set proliferation signals to 0.
										ph_[x+i][y+j] = 0;
										
										//Reset no. of proliferations remaining until death.
										rhos_[x+i][y+j] = rho;
										
										//Kill the cell if it has reached the limit.
										if (--rhos_[x][y] == 0)
										{
											cellState(x, y, DEAD);
											awakeNeighbourhood(x, y);
										}
									}
									else
									{
										//Move to the selected position
										cellState(x, y, DEAD);
										cellState(x + i, y + j, MIGRATED);
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
			cellState(x, y, DEAD);
			
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
	
	delete barrier_;
}
