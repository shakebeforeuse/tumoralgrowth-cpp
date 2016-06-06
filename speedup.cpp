#include <chrono>
#include <iostream>
#include "tumor_automaton.h"

int main(int argc, char* argv[])
{
	int size, it, tasks, step;
	
	if (argc != 5)
	{
		std::cout << "Input size" << std::endl;
		std::cin >> size;
		std::cout << "Input no. of threads to run" << std::endl;
		std::cin >> tasks;
		std::cout << "Input step between tasks" << std::endl;
		std::cin >> step;
		std::cout << "Input No. of generations" << std::endl;
		std::cin >> it;
	}
	else
	{
		size  = std::stoi(argv[1]);
		tasks = std::stoi(argv[2]);
		step  = std::stoi(argv[3]);
		it    = std::stoi(argv[4]);
	}
	
	TumorAutomaton tumor(size);
	tumor.ps  = 1;
	tumor.pp  = .8;
	tumor.pm  = .2;
	tumor.np  = 5;
	tumor.rho = 2;
	
	tumor.cellState(size/2, size/2, TumorAutomaton::ALIVE);
	
	std::cout << "Tasks\tSpeedup\tTime" << std::endl;
	
	std::chrono::time_point<std::chrono::steady_clock> tic, toc;
	
	tic = std::chrono::steady_clock::now();
	tumor.execute(it);
	toc = std::chrono::steady_clock::now();
	
	double timeSec = std::chrono::duration_cast<std::chrono::nanoseconds>(toc-tic).count();
	
	std::cout << 1 << '\t' << 1. << '\t' << timeSec * 1e-9 << std::endl;
	
	for (int i = 2; i <= tasks; i *= step)
	{
		tumor.reset();
		tumor.cellState(size/2, size/2, TumorAutomaton::ALIVE);
		
		tic = std::chrono::steady_clock::now();
		tumor.threads(i);
		tumor.execute(it);
		toc = std::chrono::steady_clock::now();
		
		double time = std::chrono::duration_cast<std::chrono::nanoseconds>(toc-tic).count();
		
		std::cout << i << '\t' << timeSec/time << '\t' << time * 1e-9 << std::endl;
	}
}
