
#include <stdexcept>
#include "photospline/I3SplineTable.h"
extern "C" {
#include "photospline/bspline.h"
}

static void
splinetable_destructor(splinetable *table) {
	if (!table) return;
	
	splinetable_free(table);
	delete table;
}

I3SplineTable::I3SplineTable(const std::string &path)
{
	table_ = boost::shared_ptr<splinetable>(new splinetable,
            splinetable_destructor);

	if (readsplinefitstable(path.c_str(), &*table_) != 0)
		throw std::runtime_error("Couldn't read spline table " + path);
}

int
I3SplineTable::Eval(double *coordinates, double *result)
{
	int centers[table_->ndim];
	
	if (tablesearchcenters(&*table_, coordinates, centers) == 0)
		*result = ndsplineeval(&*table_, coordinates, centers, 0);
	else
		return EINVAL;
	
	return 0;
}
