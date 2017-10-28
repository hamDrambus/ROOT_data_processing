#ifndef SIGNAL_OPERATIONS_H
#define SIGNAL_OPERATIONS_H

#include "GlobalDefinitions.h"

namespace SignalOperations
{

	void invert_y(std::vector<double> &x_in_out, std::vector<double> &y_in_out);
	double find_baseline_by_median(double approx, std::vector<double>&xs, std::vector<double>&ys, std::vector<peak> &peaks);
	double find_baseline_by_integral(double approx, std::vector<double>&xs, std::vector<double>&ys, std::vector<peak> &peaks);
	void find_peaks(DVECTOR &xs, DVECTOR &ys, std::vector<peak> &peaks, double base_line, double threshold, int N_trust);
	void integrate(std::vector<double>&xs, std::vector<double>&ys, std::vector<double> &y_out, double baseline = 0);
	void integrate(std::vector<double>&xs, std::vector<double>&ys, std::vector<double> &x_out, std::vector<double> &y_out,
		double left, double right, double baseline = 0);
	void integrate(DVECTOR &xs, DVECTOR &ys, DVECTOR &x_out, DVECTOR &y_out,
		DITERATOR left, DITERATOR right, double baseline = 0);
	void integrate(DVECTOR &xs, DVECTOR &ys, double &y_out,
		DITERATOR left, DITERATOR right, double baseline = 0);
	//^area is [a,b], not (a,b)
	void apply_time_limits(std::vector<double>&xs, std::vector<double>&ys, double x_left, double x_right);

	std::vector<double>::iterator find_x_iterator_by_value(std::vector<double>::iterator &x_left, std::vector<double>::iterator &x_right, double x);
	void get_max(std::vector<double>&xs, std::vector<double>&ys, std::vector<double>::iterator &x_max, double &y_max,
		int N_trust = 1);
	void find_next_peak(std::vector<double>&xs, std::vector<double>&ys, std::vector<double>::iterator &x_start,
		std::vector<double>::iterator &x_finish, double threshold, int N_trust);
	//void refine_intersection(std::vector<double>&xs, std::vector<double>&ys, std::vector<double>::iterator &x_approx, std::vector<double>::iterator minimal_x, double threshold, int N_trust);
	////for the left peak slope
	//void refine_intersection(std::vector<double>&xs, std::vector<double>&ys, std::vector<double>::iterator &x_approx, double threshold, int N_trust);
	////for the right peak slope

	void find_next_extremum(std::vector<double>&xs, std::vector<double>&ys, std::vector<double>::iterator &x_start, int N_trust);
	//TODO: theese 5 above functions rely on the 2nd order polynom properties, so it will be neater to create 2nd order polynome on the range [a,b]
	//class with such methods as getting maximum, obtaining iterators and such.
	//void find_GEM_start_time(std::vector<double>&xs, std::vector<double>&ys, std::vector<double>::iterator &x_start, int N_trust = 10);
	//N trust is the number of points over which 2nd order interpolation takes effect
	//this is quite specific and complex operation, maybe should be transferred to AnalysisManager
	void spread_peaks(double x_left, double x_right, std::vector<peak> &peaks);
	void spread_peaks(DVECTOR &xs_in, DVECTOR &ys_in, DVECTOR &xs_out, DVECTOR& ys_out);
};

#endif