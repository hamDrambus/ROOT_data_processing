#ifndef SAVITZKY_GOLAY_FILTER_H
#define SAVITZKY_GOLAY_FILTER_H

#include "GlobalParameters.h"
#ifndef _AVOID_CERN_ROOT
#include "TMatrixD.h"
#include "TVectorD.h"
#endif //_AVOID_CERN_ROOT
#include "PolynomialFit.h"
//TODO: add const S y(x)dx condition

class SavitzkyGolayFilter
{
protected:
	std::size_t _n_points;
	std::size_t _n_iterations;
	std::size_t _order;
	std::vector<double> calculate_coefs(void) const; //for equidistant points
	double pown(double val, unsigned int n) const;
public:
	SavitzkyGolayFilter(std::size_t n_points = 10, std::size_t order = 4, std::size_t n_iterations = 1);
	void setNPoints(std::size_t n);
	void setOrder(std::size_t n);
	void setNIter(std::size_t n);
	void setPars(std::size_t n_points = 10, std::size_t order = 4, std::size_t n_iterations = 1);

	std::size_t getNPoints(void) const;
	std::size_t getOrder(void) const;
	std::size_t getNIter(void) const;
	void getPars(std::size_t &n_points, std::size_t &order, std::size_t &n_iterations) const;

	void operator ()(DVECTOR &xs_in_out, DVECTOR &ys_in_out) const;
	//Assumes xs are equidistant
	void operator ()(DVECTOR &ys_in_out) const;
	bool isValid(void) const;
};

#endif
