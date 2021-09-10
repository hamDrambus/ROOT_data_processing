#include <stdio.h>

#include "GlobalParameters.h"

DITERATOR iter_add(DITERATOR& to, int what, DITERATOR& end)
{
	if (what < 0)
		return end;
	return ((int)(end - to) < what) ? end : to + what;
}

void open_output_file(std::string name, std::ofstream &str, std::ios_base::openmode _mode)
{
	ensure_file(name);
	str.open(name.c_str(), _mode);
	if (!str.is_open()){
		std::cout << "Failed to open \"" << name << "\"" << std::endl;
	}
}

void ensure_file(std::string fname)
{
	std::string folder = fname;
	while ((folder.back() != '\\') &&(folder.back()!='/') &&!folder.empty())
		folder.pop_back();
	if (!folder.empty())
		folder.pop_back();
	ensure_folder(folder);
}

void ensure_folder(std::string folder)
{
#if defined(__WIN32__)
	if (!folder.empty()) {
		DWORD ftyp = GetFileAttributesA(folder.c_str());
		if (!(ftyp & FILE_ATTRIBUTE_DIRECTORY) || ftyp == INVALID_FILE_ATTRIBUTES) {
			int code = system(("mkdir \"" + folder + "\"").c_str());
			if (code)
				std::cout << "mkdir error: " << GetLastError() << std::endl;
		}
	}
#else
	struct stat st;
	if (-1==stat(folder.c_str(), &st)) {
		int err = errno;
		switch (err) {
		case (EACCES): {
			std::cout<<"Access error"<<std::endl;
			break;
		}
		case (ENAMETOOLONG): {
			std::cout<<"Path is too long"<<std::endl;
			break;
		}
		case (ENOENT) :
		case (ENOTDIR): {
			int code = system(("mkdir -p \"" + folder + "\"").c_str());
			if (code)
				std::cout << "mkdir -p error: " << code << std::endl;
			break;
		}
		default:{
			std::cout<<"stat(\""<<folder<<"\") returned -1; errno == "<<err<<std::endl;
			break;
		}
		}
	} else {
		if (!S_ISDIR(st.st_mode)) {
			int code = system(("mkdir -p \"" + folder + "\"").c_str());
			if (code)
				std::cout << "mkdir -p error: " << code << std::endl;
		}
	}
#endif //_WIN32__
}

std::string strtoken(std::string &in, std::string break_symbs)
{
	std::string out;
	while (!in.empty()) {
		char a = in.front();
		in.erase(in.begin());
		bool break_ = false;
		for (auto h = break_symbs.begin(); h != break_symbs.end(); ++h)
			if (a == *h) {
				break_ = true;
				break;
			}
		if ((break_) && (out.empty()))
			continue;
		if (break_)
			return out;
		out.push_back(a);
	}
	return out;
}

namespace ParameterPile
{
	analysis_manifest gManifest;
	std::string this_path;
	int threads_number; //obv. must be >=1
	std::size_t max_pics_number;
	bool gnuplot_presits;
	bool quiet_mode;

	int Max_iteration_N = 1;

	int gnuplot_pad_size = 400;
	int gnuplot_max_size = 1600;
	int gnuplot_width = 1200; //default for gnuplot is 640

	bool read_accepted_events(std::string file, accepted_events<double> &info)
	{
		std::ifstream str;
		str.open(file);
		if (!str.is_open()) {
			std::cerr << "ParameterPile::read_accepted_events: Warning: Failed to open \"" << file << "\"" << std::endl;
			return false;
		}
		std::string line, word;
		int line_n = 0;
		while (!str.eof() && str.is_open()) {
			std::getline(str, line);
			++line_n;
			if (line.size() >= 2) //Ignore simple c style comment
				if ((line[0] == '/') && (line[1] == '/'))
					continue;
			try {
				word = strtoken(line, "\t ");
				if (word.empty())
					continue;
				int run = std::stoi(word);
				word = strtoken(line, "\t ");
				if (word.empty())
					continue;
				int subrun = std::stoi(word);
				word = strtoken(line, "\t ");
				double trigger = 0;
				if (!word.empty())
					trigger = std::stod(word);
				info.push(run, subrun, trigger);
			}
			catch (std::exception &e) {
				std::cerr << "ParameterPile::read_accepted_events: Error: Exception in \"" << file << "\"" << std::endl
					<< "\ton line " << line_n << std::endl;
				std::cerr << "\t" << e.what() << std::endl;
				continue;
			}
		}
		return true;
	}

	void Init_globals(void)
	{
		char path[FILENAME_MAX];
#if defined(__WIN32__)
		this_path = _getcwd(path, FILENAME_MAX);
#else
		this_path = getcwd(path, FILENAME_MAX);
#endif //__WIN32__
		if (!this_path.empty())
			if (this_path.back()!='/')
				this_path.push_back('/');

		threads_number = 8; //obv. must be >=1
		max_pics_number = 100;
		gnuplot_presits = false; //TODO: ->gnuplot_persist
		quiet_mode = true;

		//Init190404_tests(gManifest);
		//Init190404(gManifest);
		//Init180705_tests(gManifest);
        //Init180705(gManifest);
		//Init180830Xray_tests(gManifest);
		//Init180830Xray(gManifest);
		//Init191107(gManifest);
		//Init191128_tests(gManifest);
		//Init191128(gManifest);
		//Init200116_tests(gManifest);
		//Init200116(gManifest);
		//Init200213_tests(gManifest);
		//Init200213(gManifest);
		//Init170622Cd_tests(gManifest);
		//Init170622Cd(gManifest);
		//Init200910Pu_tests(gManifest);
		//Init200910Pu(gManifest);
		//Init190307Xray_tests(gManifest);
		//Init190307Xray(gManifest);
		//Init201015_tests(gManifest);
		//Init201015(gManifest);
		//Init201015Xray_tests(gManifest);
		//Init201015Xray(gManifest);
		//Init201217_tests(gManifest);
		//Init201217(gManifest);
		//Init210121_tests(gManifest);
		//Init210121(gManifest);
		//Init210128_tests(gManifest);
		//Init210128(gManifest);
		//Init210218_tests(gManifest);
		//Init210218(gManifest);
		//Init210302_tests(gManifest);
		//Init210302(gManifest);
		//Init210311_tests(gManifest);
		//Init210311(gManifest);
		//Init210316_tests(gManifest);
		//Init210316(gManifest);
		//Init210429_tests(gManifest);
		//Init210429(gManifest);
		//Init210513_tests(gManifest);
		//Init210513(gManifest);
		//Init210603_tests(gManifest);
		//Init210603(gManifest);
		//Init210902_tests(gManifest);
		Init210902(gManifest);
	}

};
