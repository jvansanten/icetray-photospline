#include <stdexcept>
#include <photospline/I3SplineTable.h>

extern "C" {
	#include <photospline/bspline.h>
}

I3SplineTable::I3SplineTable(const std::string &path)
{
	if (readsplinefitstable(path.c_str(), &table_) != 0)
		throw std::runtime_error("Couldn't read spline table " + path);
}

I3SplineTable::~I3SplineTable()
{
	splinetable_free(&table_);
}

int
I3SplineTable::Eval(double *coordinates, double *result)
{
	int centers[table_.ndim];
	
	if (tablesearchcenters(&table_, coordinates, centers) == 0)
		*result = ndsplineeval(&table_, coordinates, centers, 0);
	else
		return EINVAL;
	
	return 0;
}

