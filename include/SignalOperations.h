#ifndef SIGNAL_OPERATIONS_H
#define SIGNAL_OPERATIONS_H

#include "GlobalDefinitions.h"

namespace SignalOperations
{

	void invert_y(std::vector<double> &x_in_out, std::vector<double> &y_in_out);
	double find_baseline_by_median(double approx, std::vector<double>&xs, std::vector<double>&ys, std::vector<peak> &peaks);
	double find_baseline_by_integral(double approx, std::vector<double>&xs, std::vector<double>&ys, std::vector<peak> &peaks);
	void integrate(std::vector<double>&xs, std::vector<double>&ys, std::vector<double> &y_out, double baseline = 0);
	void integrate(std::vector<double>&xs, std::vector<double>&ys, std::vector<double> &x_out, std::vector<double> &y_out,
		double left, double right, double baseline = 0);
	void integrate(DVECTOR &xs, DVECTOR &ys, DVECTOR &x_out, DVECTOR &y_out,
		DITERATOR left, DITERATOR right, double baseline = 0);
	void integrate(DVECTOR &xs, DVECTOR &ys, double &y_out,
		DITERATOR left, DITERATOR right, double baseline = 0);
	//^area is [a,b], not (a,b)
	void apply_time_limits(std::vector<double>&xs, std::vector<double>&ys, double x_left, double x_right);

	DITERATOR find_x_iterator_by_value(DITERATOR &x_left, DITERATOR &x_right, double x);
	void get_max(DVECTOR &xs, DVECTOR &ys, DITERATOR &x_max, double &y_max, int N_trust = 1);
	void get_max(DVECTOR &xs, DVECTOR &ys, DITERATOR x_start, DITERATOR x_finish, DITERATOR &x_max, double &y_max, int N_trust = 1);
	//finds max only at [x_start, x_finish)
	void find_next_peak(DVECTOR&xs, DVECTOR&ys, DITERATOR &x_start, DITERATOR &x_finish, double threshold, int N_trust);
	//seaches peak from x_start
	void find_peaks(DVECTOR &xs, DVECTOR &ys, std::vector<peak> &peaks, double base_line, double threshold, int N_trust);

	void find_next_extremum(std::vector<double>&xs, std::vector<double>&ys, std::vector<double>::iterator &x_start, int N_trust);
	//N trust is the number of points over which 2nd order interpolation takes effect
	//done: TODO: theese 5 above functions rely on the 2nd order polynom properties, so it will be neater to create 2nd order polynome on the range [a,b]
	//class with such methods as getting maximum, obtaining iterators and such.
	void spread_peaks(double x_left, double x_right, std::vector<peak> &peaks, DVECTOR &xs_out, DVECTOR& ys_out);
	void peaks_to_yx(double x_left, double x_right, std::vector<peak> &peaks, DVECTOR &xs_out, DVECTOR& ys_out);
	void spread_peaks(DVECTOR &xs_in, DVECTOR &ys_in, DVECTOR &xs_out, DVECTOR& ys_out);

	void substract_baseline(DVECTOR &ys_in, double base_line);
	void substract_baseline(DVECTOR &ys_in, DVECTOR &base_ys);
};

#endif