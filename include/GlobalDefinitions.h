#ifndef GLOBAL_DEFINITIONS_H
#define GLOBAL_DEFINITIONS_H

#include <iostream>
#include <fstream>
#include <vector>
#include <deque>
#include <string.h>
#include <sstream>

#include <TThread.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TVector.h>
#include <TF1.h>
#include <TMath.h>
#include <Math/Point2D.h>
#include <windows.h>

#undef max
#undef min

#define DATA_PREFIX std::string("../../../Data/")
#define DATA_NAME_FORMAT "^run_\d+__ch_\d+\.dat$"
#define DATA_EXPERIMENT_FORMAT "^event_x-ray_\d{1,2}_.*$"

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
#define OUTPUT_DIR "results\\"
#define OUTPUT_GEMS "GEMs.txt"
#define OUTPUT_MPPCS "MPPC_"
#define OUTPUT_MPPCS_PICS "MPPCs_v2\\MPPCs_"
#define _TEMP_CODE
#define _HOTFIX_DECREASE_MPPC_MEMORY_USAGE
#define _HOTFIX_CLEAR_MEMORY

#define STD_CONT std::deque

//#define _USE_DEQUE
#ifdef _USE_DEQUE
#define DVECTOR std::deque<double>
#else
#define DVECTOR std::vector<double>
#endif

#define DITERATOR DVECTOR::iterator
#define D_REV_ITERATOR DVECTOR::reverse_iterator

//#define _PROCESS_GEMS

void open_output_file(std::string name, std::ofstream &str);

struct peak //TODO: add amplitude
{
	double left;
	double right;
	double S; //Area
	double A; //Amplitude (from baseline)
};

#endif