
#ifndef PHOTOSPLINE_I3SPLINETABLE_H_INCLUDED
#define PHOTOSPLINE_I3SPLINETABLE_H_INCLUDED

#include <string>
#include <photospline/splinetable.h>

class I3SplineTable {
public:
	/**
	 * @param[in] path Path to a FITS file
	 */ 
	I3SplineTable(const std::string &path);
	virtual ~I3SplineTable();

	/** Evaluate the spline surface
	 * 
	 * @param[in]       x N-dimensonal coordinates at which to evaluate
	 * @param[out] result Value of spline surface at coordinates
	 * @returns 0 on success, non-zero otherwise
	 */
	int Eval(double *x, double *result) const;

	/** Get the number of dimensions */
	unsigned GetNDim() const { return table_.ndim; };
	/** Get the extent of full support in dimension dim */
	std::pair<double, double> GetExtents(int dim) const;
	
	/** Retrieve a value stored in the FITS header */
	double GetField(const std::string &key) const;
private:
	I3SplineTable(const I3SplineTable&);
	
	struct splinetable table_;
	double bias_;
};

#endif
