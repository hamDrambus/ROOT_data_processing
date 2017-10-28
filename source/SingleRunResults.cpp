#include "SingleRunResults.h"

SingleRunResults::SingleRunResults(SingleRunData *_of_what)
{
	of_what = _of_what;
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