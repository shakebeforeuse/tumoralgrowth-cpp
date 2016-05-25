#include <chrono>
#include <iostream>
#include "tumor_automaton.h"

int main(int argc, char* argv[])
{
	int size, it, tasks;
	
	if (argc != 4)
	{
		std::cout << "Input size" << std::endl;
		std::cin >> size;
		std::cout << "Input No. of generations" << std::endl;
		std::cin >> it;
		std::cout << "Input no. of threads to run" << std::endl;
		std::cin >> tasks;
	}
	else
	{
		size  = std::stoi(argv[1]);
		it    = std::stoi(argv[2]);
		tasks = std::stoi(argv[3]);
	}
	
	TumorAutomaton tumor(size);
	tumor.ps = 1;
	tumor.pp = 1;
	tumor.cellState(size/2, size/2, TumorAutomaton::ALIVE);
	
	std::chrono::time_point<std::chrono::steady_clock> tic, toc;
	
	tic = std::chrono::steady_clock::now();
	tumor.threads(tasks);
	tumor.execute(it);
	toc = std::chrono::steady_clock::now();
	
	double elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(toc-tic).count();
	
	std::cout << elapsed * 1e-9 << std::endl;
	
	for (int i = 0; i < size; ++i)
	{
		for (int j = 0; j < size; ++j)
			switch (tumor.cellState(i, j))
			{
				case 0:
					std::cerr << " ";
					break;
				
				case 1:
					std::cerr << ".";
					break;
					
				case 2:
					std::cerr << "o";
					break;
				
				case 3:
					std::cerr << "O";
					break;
				
				case 4:
					std::cerr << "M";
					break;
				
				default:
					std::cerr << "#";
					break;
			}
		
		std::cerr << std::endl;
	}
		
}
