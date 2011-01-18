
extern "C" {
	#include "photospline/splinetable.h"
	#include "photospline/bspline.h"
}

#include <I3Test.h>
#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;

struct TableSet {
	fs::path abs, prob;
};

TableSet
get_splinetables()
{
	ENSURE(getenv("I3_SRC") != NULL,
	    "I3_SRC must be defined in the parent shell.");

	const std::string I3_SRC(getenv("I3_SRC"));
	
	fs::path abs_table(I3_SRC + "/photospline/private/test/ems_z0_a0.pt.abs.fits");
	fs::path prob_table(I3_SRC + "/photospline/private/test/ems_z0_a0.pt.prob.fits");
		
	ENSURE(fs::exists(abs_table), "Amplitude table exists.");
	ENSURE(fs::exists(prob_table), "Quantile table exists.");
	
	TableSet tabset;
	tabset.abs = abs_table;
	tabset.prob = prob_table;
	
	return tabset;
}

TEST_GROUP(BoundaryIssues);

void splinetable_destructor(struct splinetable *table) {
	if (!table) return;
	
	splinetable_free(table);
	delete table;
	
}

TEST(QuantileContinuity)
{
	TableSet tables = get_splinetables();
	boost::shared_ptr<struct splinetable> table(new struct splinetable, splinetable_destructor);
	const unsigned time_dim = 3;
	
	ENSURE(readsplinefitstable(tables.prob.string().c_str(), table.get()) == 0, "Table can be read.");
	
	unsigned i;
	unsigned ndim = table->ndim;
	double tablecoords[ndim], base, nudge;
	int centers[ndim];
	
	for (i = 0; i < ndim; i++) {
		double low, high;
		/* Don't trust the extents for the time dimension; they lie. */
		if (i == time_dim) {
			low = table->knots[i][table->order[i]];
			high = table->knots[i][table->naxes[i]];
			tablecoords[i] = low - std::numeric_limits<double>::epsilon();
		} else {
			low = table->extents[i][0];
			high = table->extents[i][1];
			tablecoords[i] = (low+high)/2.0;
		}
	}
	
	ENSURE(tablecoords[time_dim] > 0.0, "t=0 is not fully supported.");
	
	ENSURE(tablesearchcenters(table.get(), tablecoords, centers) == 0,
	    "tablesearchcenters() succeeds for t>0.");
	
	base = ndsplineeval(table.get(), tablecoords, centers, 0);
	
	ENSURE(base >= 0, "Base quantile is positive.");
	
	tablecoords[time_dim] += std::numeric_limits<double>::epsilon();
	
	ENSURE(tablesearchcenters(table.get(), tablecoords, centers) == 0,
	    "tablesearchcenters() still succeeds.");
	
	nudge = ndsplineeval(table.get(), tablecoords, centers, 0);
	
	ENSURE(nudge >= 0, "Quantile at edge+epsilon is also positive.");
	
	ENSURE_DISTANCE(base, nudge, std::numeric_limits<float>::epsilon(),
	    "Time quantile is continuous at the edge of support "
	    "(slope is less than FLOAT_EPSILON/DBL_EPSILON)");
}