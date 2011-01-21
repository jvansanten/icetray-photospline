
extern "C" {
	#include "photospline/splinetable.h"
	#include "photospline/bspline.h"
}

#include <I3Test.h>
#include <boost/filesystem.hpp>
#include <sys/time.h>

namespace fs = boost::filesystem;

struct TableSet {
	fs::path abs, prob;
};

static void
splinetable_destructor(struct splinetable *table) {
	if (!table) return;
	
	splinetable_free(table);
	delete table;
	
}

static TableSet
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

static boost::shared_ptr<struct splinetable>
load_splinetable(fs::path &fname)
{
	boost::shared_ptr<struct splinetable> table(new struct splinetable, splinetable_destructor);
	ENSURE(readsplinefitstable(fname.string().c_str(), table.get()) == 0, "Table can be read.");
	
	return table;
}


TEST_GROUP(BoundaryIssues);

TEST(QuantileContinuity)
{
	const unsigned time_dim = 3;
	TableSet tables = get_splinetables();
	boost::shared_ptr<struct splinetable> table = load_splinetable(tables.prob);
	
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

/*
 * bsplvb_simple() can be made to return sensical values anywhere
 * in the knot field.
 */

TEST(bsplvb_simple_vs_bspline)
{
	unsigned i;
	const unsigned n_knots = 10;
	const int order = 2;
	double x, *knots;
	int center, offset;
	float localbasis_bsplvb[order+1], localbasis_bspline[order+1];
	struct timeval thetime;
	std::vector<double> knotvec;
	
	/* Seed our crappy little PSRNG. */
	gettimeofday(&thetime, NULL);
	srand(thetime.tv_sec);
	/* Generate a random knot field. */
	for (i = 0; i < n_knots; i++)
		knotvec.push_back(double(rand())/double(RAND_MAX));
	std::sort(knotvec.begin(), knotvec.end());
	knots = &knotvec.front();
	
	/* Before the first fully-supported knot. */
	for (i = 0; i < order+1; i++) {
		x = (knots[i]+knots[i+1])/2.0; /* Less than fully-supported */
		center = order; /* First fully-supported spline. */
	
		ENSURE(int(i) <= center);
	
		bsplvb_simple(knots, n_knots, x, center /* left */,
		    order+1 /* jhigh */, localbasis_bsplvb /* biatx */);
	
		for (offset = -order; offset <= 0; offset++) {
			ENSURE(offset+order >= 0);
			localbasis_bspline[offset+order] =
			    bspline(knots, x, center + offset, order);
		}
	
		for (offset = 0; offset < order+1; offset++) {
			ENSURE_DISTANCE(localbasis_bspline[offset], localbasis_bsplvb[offset],
			std::numeric_limits<double>::epsilon());
		}
	}
	
	/* Within the support. */
	for (i = order+1; i < n_knots-order-1; i++) {
		x = (knots[i]+knots[i+1])/2.0;
		center = i;
		
		bsplvb_simple(knots, n_knots, x, center /* left */,
		    order+1 /* jhigh */, localbasis_bsplvb /* biatx */);
	
		for (offset = -order; offset <= 0; offset++) {
			ENSURE(offset+order >= 0);
			localbasis_bspline[offset+order] =
			    bspline(knots, x, center + offset, order);
		}
	
		for (offset = 0; offset < order+1; offset++) {
			ENSURE_DISTANCE(localbasis_bspline[offset], localbasis_bsplvb[offset],
			std::numeric_limits<double>::epsilon());
		}
	}
	
	/* After the last first fully-supported knot. */
	for (i = n_knots-order-1; i < n_knots-2; i++) {
		x = (knots[i]+knots[i+1])/2.0; /* Less than fully-supported */
		center = n_knots-order-2; /* Last fully-supported spline. */
	
		ENSURE(int(i) >= center);
	
		bsplvb_simple(knots, n_knots, x, center /* left */,
		    order+1 /* jhigh */, localbasis_bsplvb /* biatx */);
	
		for (offset = -order; offset <= 0; offset++) {
			ENSURE(offset+order >= 0);
			localbasis_bspline[offset+order] =
			    bspline(knots, x, center + offset, order);
		}
	
		for (offset = 0; offset < order+1; offset++) {
			ENSURE_DISTANCE(localbasis_bspline[offset], localbasis_bsplvb[offset],
			std::numeric_limits<double>::epsilon());
		}
	}
}

/*
 * bspline_deriv_nonzero() can be made to return sensical values anywhere
 * in the knot field.
 */

TEST(bspline_deriv_nonzero_vs_bspline_deriv)
{
	unsigned i;
	const unsigned n_knots = 10;
	const int order = 2;
	double x, *knots;
	int center, offset;
	float localbasis_bsplvb[order+1], localbasis_bspline[order+1];
	struct timeval thetime;
	std::vector<double> knotvec;
	/* This calculation is less stable. */
	double tol = 10*std::numeric_limits<float>::epsilon();
	
	/* Seed our crappy little PSRNG. */
	gettimeofday(&thetime, NULL);
	srand(thetime.tv_sec);
	/* Generate a random knot field. */
	for (i = 0; i < n_knots; i++)
		knotvec.push_back(double(rand())/double(RAND_MAX));
	std::sort(knotvec.begin(), knotvec.end());
	knots = &knotvec.front();
	
	/* Before the first fully-supported knot. */
	for (i = 0; i < order+1; i++) {
		x = (knots[i]+knots[i+1])/2.0; /* Less than fully-supported */
		center = order; /* First fully-supported spline. */
	
		ENSURE(int(i) <= center);
	
		bspline_deriv_nonzero(knots, n_knots, x, center /* left */,
		    order /* jhigh */, localbasis_bsplvb /* biatx */);
	
		for (offset = -order; offset <= 0; offset++) {
			ENSURE(offset+order >= 0);
			localbasis_bspline[offset+order] =
			    bspline_deriv(knots, x, center + offset, order);
		}
	
		for (offset = 0; offset < order+1; offset++) {
			ENSURE_DISTANCE(localbasis_bspline[offset], localbasis_bsplvb[offset], tol);
		}
	}
	
	/* Within the support. */
	for (i = order+1; i < n_knots-order-1; i++) {
		x = (knots[i]+knots[i+1])/2.0;
		center = i;
		
		bspline_deriv_nonzero(knots, n_knots, x, center /* left */,
		    order, localbasis_bsplvb /* biatx */);
	
		for (offset = -order; offset <= 0; offset++) {
			ENSURE(offset+order >= 0);
			localbasis_bspline[offset+order] =
			    bspline_deriv(knots, x, center + offset, order);
		}
	
		for (offset = 0; offset < order+1; offset++) {
			ENSURE_DISTANCE(localbasis_bspline[offset], localbasis_bsplvb[offset], tol);
		}
	}

	/* After the last first fully-supported knot. */
	for (i = n_knots-order-1; i < n_knots-2; i++) {
		x = (knots[i]+knots[i+1])/2.0; /* Less than fully-supported */
		center = n_knots-order-2; /* Last fully-supported spline. */
	
		ENSURE(int(i) >= center);
	
		bspline_deriv_nonzero(knots, n_knots, x, center /* left */,
		    order /* jhigh */, localbasis_bsplvb /* biatx */);
	
		for (offset = -order; offset <= 0; offset++) {
			ENSURE(offset+order >= 0);
			localbasis_bspline[offset+order] =
			    bspline_deriv(knots, x, center + offset, order);
		}
	
		for (offset = 0; offset < order+1; offset++) {
			ENSURE_DISTANCE(localbasis_bspline[offset], localbasis_bsplvb[offset], tol);
		}
	}
}