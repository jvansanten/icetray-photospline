
#ifndef PHOTOSPLINE_I3SPLINETABLE_H_INCLUDED
#define PHOTOSPLINE_I3SPLINETABLE_H_INCLUDED

#include <string>

extern "C" {
	#include <photospline/splinetable.h>
}

class I3SplineTable {
public:
	I3SplineTable(const std::string &path);
	virtual ~I3SplineTable();

	int Eval(double *x, double *result) const;

	unsigned GetNDim() const { return table_.ndim; };
	std::pair<double, double> GetExtents(int dim) const;
private:
	struct splinetable table_;
	double bias_;
};

#endif
