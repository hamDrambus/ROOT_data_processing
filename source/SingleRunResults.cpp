#include "SingleRunResults.h"

SingleRunResults::SingleRunResults(SingleRunData *_of_what)
{
	curr_area = _of_what ? _of_what->getArea():curr_area;
	setValid(false);
	_current_status = Status::Empty;
	PMT3_n_peaks = 0;
	PMT3_summed_peaks_area = 0;
}

SingleRunResults::Status SingleRunResults::getStatus(void) const
{	return _current_status;}

bool SingleRunResults::isValid(void) const
{	return is_valid;}
void SingleRunResults::setValid(bool val)
{	is_valid = val;}

std::size_t SingleRunResults::real_size(void)
{
	std::size_t results_size = sizeof(*this);
	results_size += xs_GEM.capacity()*sizeof(double);
	results_size += ys_GEM.capacity()*sizeof(double);

	results_size += mppc_peaks.size()*sizeof(STD_CONT<peak>);
	for (int ch=0, _end_ch=mppc_peaks.size(); ch!=_end_ch; ++ch)
		results_size += mppc_peaks[ch].size()*sizeof(peak);
	results_size += pmt_peaks.size()*sizeof(STD_CONT<peak>);
	for (int ch=0, _end_ch=pmt_peaks.size(); ch!=_end_ch; ++ch)
		results_size += pmt_peaks[ch].size()*sizeof(peak);

	results_size += PMT_S2_integral.capacity()*sizeof(DVECTOR);

	results_size += mppc_baseline_xs.size()*sizeof(DVECTOR);
	for (int ch=0, _end_ch=mppc_baseline_xs.size(); ch!=_end_ch; ++ch)
		results_size += mppc_baseline_xs[ch].size()*sizeof(peak);

	results_size += mppc_baseline_ys.size()*sizeof(DVECTOR);
	for (int ch=0, _end_ch=mppc_baseline_ys.size(); ch!=_end_ch; ++ch)
		results_size += mppc_baseline_ys[ch].size()*sizeof(peak);

	results_size += mppc_S2_peaks_area.capacity()*sizeof(double);
	results_size += mppc_S2_start_t.capacity()*sizeof(double);
	results_size += mppc_S2_finish_t.capacity()*sizeof(double);
	results_size += mppc_double_I.capacity()*sizeof(double);

	results_size += mppc_channels.size()*sizeof(int);
	results_size += pmt_channels.size()*sizeof(int);

	results_size += xs_PMT3.capacity()*sizeof(double);
	results_size += ys_PMT3.capacity()*sizeof(double);
	results_size += xs_PMT1.capacity()*sizeof(double);
	results_size += ys_PMT1.capacity()*sizeof(double);
	return results_size;
}
