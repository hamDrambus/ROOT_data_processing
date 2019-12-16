#include <ctime>
#include "GlobalParameters.h"
#include "MTAnalysisManager.h"

//13.12.2019 TODO: so much rework is required
//	0) add description of program structure as well as algorithms
//	1) ParameterPile--->gSettings (renaming)
//	2) make area_vector to be thread-save
//	3) support writing data in ROOT trees, not only in my custom files (too cumbersome hierarchy)
//	*4) Qt (Python?) interactive signal displaying and on-the run selection of processing parameters (like oscilloscope)
//	5) Improve readability of SignalOperations
//	!!!6) Algorithm for separation of slightly merged peaks f(baseline, threshold, running_threshold = fraction of peak maximum)
//	7) add Vlad's algorithm for peak search

int main(int argc, char *argv[])
{
	ParameterPile::Init_globals();
	std::cout << "Working directory: " << ParameterPile::this_path << std::endl;
	
#ifdef _USE_TIME_STATISTICS
	auto start = std::chrono::system_clock::now();
	std::cout << "Starting processing data on ";
	std::time_t start_time = std::chrono::system_clock::to_time_t(start);
	std::cout << std::ctime(&start_time) << std::endl;
#endif
	if (1 >= ParameterPile::threads_number) {
		AnalysisManager man(ParameterPile::gManifest);
		man.processAllExperiments();
	} else {
		MTAnalysisManager man(ParameterPile::gManifest);
		man.processAllExperiments();
	}
	std::cout<<std::endl<<"========================================="<<std::endl << "Finished" << std::endl;
#ifdef _USE_TIME_STATISTICS
	std::cout << "Started processing data on ";
	std::cout << std::ctime(&start_time) << std::endl;
	auto end = std::chrono::system_clock::now();
	auto diff = end - start;
	std::cout << "Finished processing on ";
	start_time = std::chrono::system_clock::to_time_t(end);
	std::cout << std::ctime(&start_time) << std::endl;
	std::cout << "Elapsed time  = \t" << std::chrono::duration_cast<std::chrono::seconds>(diff).count() << " s."<< std::endl;
#endif
	return 0;
}
