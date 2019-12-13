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
#include <chrono>
#if defined(__WIN32__)
#include <sehmap.h>
#include <direct.h>
#include <Windows4Root.h>
#else
#include <dirent.h>
#include <errno.h>
#include <sys/stat.h>
#endif
#include <limits>
#include <float.h>

#include <TThread.h>
#include <TSpectrum.h>
#include <TApplication.h>
#include <TMath.h>

#undef max
#undef min

#define _USE_TIME_STATISTICS

#define STD_CONT std::deque

//#define _USE_DEQUE
#ifdef _USE_DEQUE
#define DVECTOR std::deque<double>
#else
#define DVECTOR std::vector<double>
#endif

#define DITERATOR DVECTOR::iterator
#define D_REV_ITERATOR DVECTOR::reverse_iterator

#define GET_MACRO(_1,_2, NAME,...) NAME
#define INVOKE_GNUPLOT(...) GET_MACRO(__VA_ARGS__, INVOKE_GNUPLOT2, INVOKE_GNUPLOT1)(__VA_ARGS__)

#if defined(__WIN32__)
#define INVOKE_GNUPLOT1(a) system(("start \"\" \"%GNUPLOT%\\gnuplot.exe\" --persist \"" + a + "\"").c_str())
#else
#define INVOKE_GNUPLOT1(a) system(("gnome-terminal -- bash -c \"gnuplot \"" + a +"\"\"").c_str());
#define INVOKE_GNUPLOT2(a, b) system(("gnome-terminal -- bash -c \"cd \""+b+"\"; gnuplot \"" + a +"\"\"").c_str());
#endif //__WIN32__


DITERATOR iter_add(DITERATOR& to, int what, DITERATOR& end);
void open_output_file(std::string name, std::ofstream &str, std::ios_base::openmode _mode = std::ios_base::trunc);
void ensure_file(std::string fname); //makes sure file can be created later on
void ensure_folder(std::string folder);
std::string strtoken(std::string &in, std::string break_symbs);
bool isSameChannels(const STD_CONT<int>& a, const STD_CONT<int>& b);
int getIndex(const STD_CONT<int>& channels, int ch);

class peak
{
public:
	double left;
	double right;
	double S; //Area
	double A; //Amplitude (from baseline)
	double t;
	peak();
};

#endif
