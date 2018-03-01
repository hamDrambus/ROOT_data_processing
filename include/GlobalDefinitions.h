#ifndef GLOBAL_DEFINITIONS_H
#define GLOBAL_DEFINITIONS_H

//#define __WIN32__

#include <iostream>
#include <fstream>
#include <vector>
#include <deque>
#include <string.h>
#include <sstream>
#include <functional>
#if defined(__WIN32__)
#include <direct.h>
#else
#include <dirent.h>
#include <errno.h>
#include <sys/stat.h>
#endif


#include <TThread.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TVector.h>
#include <TF1.h>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <Math/Point2D.h>

#undef max
#undef min

#define DATA_PREFIX std::string("../../Data/170622/")
#define DATA_NAME_FORMAT "^run_\d+__ch_\d+\.dat$"
#define DATA_EXPERIMENT_FORMAT "^X-ray_\d{1,2}_.*$"

#define RUNS_CACHE_PATH "output/processed_runs.txt"
#define ALL_RUNS_PATH "output/processed_experiments.txt"
#define ALL_EXPERIMENTS_PATH "output/processed_summary.txt"

#define DATA_TIME_CONSTANT 1.6e-2
//^in microseconds
#define DATA_VOLTAGE_CHANNELS 4095
#define DATA_VOLTAGE_AMPLITUDE 2.0
//in volts
#define DATA_VOLTAGE_OF_ZERO_CHANNEL (-1.0)
//in volts
#define OUTPUT_DIR std::string("../../Data/170622/results/")
//GEM_v1 - finding baseline for every event
//GEM_v2 - finding baseline for averaged signal
#define _PROCESS_GEMS
#define GEM_V2_
#undef GEM_V1_
#ifdef GEM_V1_
#define OUTPUT_GEMS "GEM_v1"
#endif
#ifdef GEM_V2_
#define OUTPUT_GEMS "GEM_v2"
#endif

#define OUTPUT_PMTS "PMT_v1/PMT_"
#define OUTPUT_MPPCS "MPPC_"
#define OUTPUT_MPPCS_PICS "MPPCs_v3/MPPCs_"
#define _TEMP_CODE
#define _HOTFIX_DECREASE_MPPC_MEMORY_USAGE
#define _HOTFIX_CLEAR_MEMORY
#define _NO_PMT_SELECTION
#define _USE_TIME_STATISTICS
//#define _DRAW_CLUSTER_FINDING

#define STD_CONT std::deque

//#define _USE_DEQUE
#ifdef _USE_DEQUE
#define DVECTOR std::deque<double>
#else
#define DVECTOR std::vector<double>
#endif

#define DITERATOR DVECTOR::iterator
#define D_REV_ITERATOR DVECTOR::reverse_iterator

#if defined(__WIN32__)
#define INVOKE_GNUPLOT(a) system(("start \"\" \"%GNUPLOT%\\gnuplot.exe\" --persist \"" + a + "\"").c_str())
#else
#define INVOKE_GNUPLOT(a) system(("konsole -e gnuplot \"" + a +"\"").c_str());
#endif //__WIN32__


DITERATOR iter_add(DITERATOR& to, int what, DITERATOR& end);
void open_output_file(std::string name, std::ofstream &str, std::ios_base::openmode _mode = std::ios_base::trunc);

class peak
{
public:
	double left;
	double right;
	double S; //Area
	double A; //Amplitude (from baseline)
	peak();
};

#endif
