#ifndef SIGNAL_OPERATIONS_H
#define SIGNAL_OPERATIONS_H

#include "GlobalParameters.h"

namespace SignalOperations
{

	void invert_y(DVECTOR &x_in_out, DVECTOR &y_in_out);
	double find_baseline_by_median(double approx, DVECTOR &xs, DVECTOR &ys, STD_CONT<peak> &peaks);
	double find_baseline_by_integral(double approx, DVECTOR &xs, DVECTOR &ys);
	double find_baseline_by_integral(double approx, DVECTOR &xs, DVECTOR &ys, STD_CONT<peak> &peaks);
	void find_baseline_by_ROOT (DVECTOR &xs, DVECTOR &ys, DVECTOR &ys_out);

	void integrate(DVECTOR &xs, DVECTOR &ys, DVECTOR &y_out, double baseline = 0);
	void integrate(DVECTOR &xs, DVECTOR &ys, DVECTOR &x_out, DVECTOR &y_out,
		double left, double right, double baseline = 0);
	void integrate(DVECTOR &xs, DVECTOR &ys, DVECTOR &x_out, DVECTOR &y_out,
		DITERATOR left, DITERATOR right, double baseline = 0);
	void integrate(DVECTOR &xs, DVECTOR &ys, double &y_out,
		DITERATOR left, DITERATOR right, double baseline = 0);
	//^area is [a,b], not (a,b)
	void apply_time_limits(DVECTOR &xs, DVECTOR &ys, double x_left, double x_right);

	DITERATOR find_x_iterator_by_value(DITERATOR &x_left, DITERATOR &x_right, double x);
	void get_max(DVECTOR &xs, DVECTOR &ys, DITERATOR &x_max, double &y_max, int N_trust = 1);
	void get_max(DVECTOR &xs, DVECTOR &ys, DITERATOR x_start, DITERATOR x_finish, DITERATOR &x_max, double &y_max, int N_trust = 1);
	//finds max only at [x_start, x_finish)
	void find_next_peak(DVECTOR&xs, DVECTOR&ys, DITERATOR &x_start, DITERATOR &x_finish, double threshold, int N_trust);
	//seaches peak from x_start
	void find_peaks(DVECTOR &xs, DVECTOR &ys, STD_CONT<peak> &peaks, double base_line, double threshold, int N_trust);
	//^these two aren't fully right functions, but ok for selection of runs by PMT signal i.e. coarse operations
	void find_next_peak_fine (DVECTOR&xs, DVECTOR&ys, DITERATOR &x_start, DITERATOR &x_finish, double &Amp,
		double thresh_finder, double thresh_edges, int N_trust);
	//seaches peak from x_start, in difference to find_next peak this one first finds peak by threshold_finder and then finds its edges (wider 
	//than intersection of threshold and signal) using thresh_edges.
	void find_peaks_fine(DVECTOR &xs, DVECTOR &ys, STD_CONT<peak> &peaks, double base_line, double threshold, double threshold_edges, int N_trust);

	void find_next_extremum(DVECTOR &xs, DVECTOR &ys, DITERATOR &x_start, int N_trust);
	//N trust is the number of points over which 2nd order interpolation takes effect
	//done: TODO: theese 5 above functions rely on the 2nd order polynom properties, so it will be neater to create 2nd order polynome on the range [a,b]
	//class with such methods as getting maximum, obtaining iterators and such.
	void spread_peaks(double x_left, double x_right, STD_CONT<peak> &peaks, DVECTOR &xs_out, DVECTOR& ys_out);
	//in comparisson to ^ has more smooth result
	void spread_peaks_v2(double x_left, double x_right, STD_CONT<peak> &peaks, DVECTOR &xs_out, DVECTOR& ys_out, double min_dx);
	void peaks_to_yx(double x_left, double x_right, STD_CONT<peak> &peaks, DVECTOR &xs_out, DVECTOR& ys_out);
	void spread_peaks(DVECTOR &xs_in, DVECTOR &ys_in, DVECTOR &xs_out, DVECTOR& ys_out);

	void exclude_peaks(DVECTOR &xs_in, DVECTOR &ys_in, STD_CONT<peak> &peaks);

	void substract_baseline(DVECTOR &ys_in, double base_line);
	void substract_baseline(DVECTOR &ys_in, DVECTOR &base_ys);
};

#endif