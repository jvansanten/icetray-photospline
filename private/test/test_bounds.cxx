
extern "C" {
	#include "photospline/splinetable.h"
	#include "photospline/bspline.h"
}

#include <I3Test.h>
#include <boost/filesystem.hpp>
#include <sys/time.h>
#include <limits>

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

/*
 * Check that analytic convolution works and preserves the monotonicity
 * of the arrival-time CDF.
 */
TEST(Convolution)
{
	TableSet tables = get_splinetables();
	boost::shared_ptr<struct splinetable> table = load_splinetable(tables.prob);
	const unsigned time_dim = 3;
	double sigma = 10.0;
	int n_knots = 3;
	double knots[3] = {-2*sigma, 0, 2*sigma};
	unsigned i;
	unsigned ndim = table->ndim;
	double tablecoords[ndim], base, nudge, q0_raw, q0_conv, q1_raw, q1_conv;
	const double eps = std::numeric_limits<double>::epsilon();
	const double tmax = table->extents[time_dim][1];
	int centers[ndim];
	
	for (i = 0; i < ndim; i++) {
		double low, high;
		low = table->extents[i][0];
		high = table->extents[i][1];
		tablecoords[i] = (low+high)/2.0;
	}
	
	tablecoords[time_dim] = table->knots[time_dim][0]+eps;
	ENSURE_EQUAL(tablesearchcenters(table.get(), tablecoords, centers), 0);
	q0_raw = ndsplineeval(table.get(), tablecoords, centers, 0);
	tablecoords[time_dim] = table->extents[time_dim][1];
	ENSURE_EQUAL(tablesearchcenters(table.get(), tablecoords, centers), 0);
	q1_raw = ndsplineeval(table.get(), tablecoords, centers, 0);
	
	ENSURE_DISTANCE(q1_raw-q0_raw, 1, 1e-3,
	    "Arrival-time CDF is normalized to within 0.1%%.");
	
	base = q0_raw;
	for (i = 1; i+table->order[time_dim]+1 < table->nknots[time_dim]; i++) {
		tablecoords[time_dim] = table->knots[time_dim][i];
		ENSURE_EQUAL(tablesearchcenters(table.get(), tablecoords, centers), 0);
		nudge = ndsplineeval(table.get(), tablecoords, centers, 0);
		ENSURE(nudge - base >= 0,
		    "Arrival-time CDF is monotonic.");
		base = nudge;
	}
	
	int err = splinetable_convolve(table.get(), time_dim, knots, n_knots);
	ENSURE_EQUAL(err, 0, "Convolution succeeds.");
	
	tablecoords[time_dim] = table->knots[time_dim][0]*(1-eps);
	ENSURE_EQUAL(tablesearchcenters(table.get(), tablecoords, centers), 0);
	q0_conv = ndsplineeval(table.get(), tablecoords, centers, 0);

	base = q0_conv;
	for (i = 1; i < table->nknots[time_dim]; i++) {
		tablecoords[time_dim] = table->knots[time_dim][i];
		if (tablecoords[time_dim] > table->extents[time_dim][1])
			break;
		ENSURE_EQUAL(tablesearchcenters(table.get(), tablecoords, centers), 0);
		nudge = ndsplineeval(table.get(), tablecoords, centers, 0);
		ENSURE(nudge - base > -1e-3,
		    "Arrival-time CDF remains monotonic.");
		base = nudge;
	}
	
	tablecoords[time_dim] = table->extents[time_dim][1];
	ENSURE_EQUAL(tablesearchcenters(table.get(), tablecoords, centers), 0);
	q1_conv = ndsplineeval(table.get(), tablecoords, centers, 0);
	
	ENSURE(table->extents[time_dim][1] < tmax,
	    "t_max is slightly smaller after convolution.");
	ENSURE(q1_conv - q1_raw > -10*std::numeric_limits<float>::epsilon(),
	    "Maximum quantile (at the end of support) is not diminished.");
	
	ENSURE_DISTANCE(q1_raw-q0_raw, q1_conv-q0_conv, 10*std::numeric_limits<float>::epsilon(),
	    "Arrival-time CDF remains normalized after convolution.");
}

/*
 * The quantiles in the time dimension are continous. This can only be true
 * if the spline evaluation code is capable of handling points near the
 * edges of the knot fields that are supported by fewer than (order+1)
 * splines.
 */
TEST(QuantileContinuity)
{
	const unsigned time_dim = 3;
	TableSet tables = get_splinetables();
	boost::shared_ptr<struct splinetable> table = load_splinetable(tables.prob);
	
	unsigned i;
	unsigned ndim = table->ndim;
	double tablecoords[ndim], base, nudge;
	const double eps = std::numeric_limits<double>::epsilon();
	int centers[ndim];
	
	for (i = 0; i < ndim; i++) {
		double low, high;
		low = table->extents[i][0];
		high = table->extents[i][1];
		tablecoords[i] = (low+high)/2.0;
	}
	
	/* Check the transition into the knot field from the left. */
	tablecoords[time_dim] = table->knots[time_dim][0];
	ENSURE(tablesearchcenters(table.get(), tablecoords, centers) != 0,
	    "tablesearchcenters() fails right at the left edge of the knot field.");
	
	tablecoords[time_dim] += eps;
	
	ENSURE_EQUAL(tablesearchcenters(table.get(), tablecoords, centers), 0,
	    "tablesearchcenters() succeeds just inside the left edge of the knot field.");
	ENSURE_EQUAL(centers[time_dim], table->order[time_dim],
	    "centers[time_dim] holds the first fully-supported knot index.");
	
	base = ndsplineeval(table.get(), tablecoords, centers, 0);
	
	ENSURE(base >= 0, "Base quantile is positive.");
	ENSURE_DISTANCE(0, base, std::numeric_limits<float>::epsilon(),
	    "Time quantile is continuous "
	    "(slope is less than FLOAT_EPSILON/DBL_EPSILON)");
		
	/*
	 * Now, step over the intertior knots, checking for continuity
	 * at every knot crossing.
	 */
	for (i = 1; i+1 < table->nknots[time_dim]; i++) {
		tablecoords[time_dim] = table->knots[time_dim][i]*(1-eps);
		ENSURE_EQUAL(tablesearchcenters(table.get(), tablecoords, centers), 0,
		    "tablesearchcenters() succeeds inside the knot field.");
		
		base = ndsplineeval(table.get(), tablecoords, centers, 0);
		ENSURE(base >= 0);
		
		tablecoords[time_dim] = table->knots[time_dim][i];
		ENSURE_EQUAL(tablesearchcenters(table.get(), tablecoords, centers), 0,
		    "tablesearchcenters() succeeds inside the knot field.");
		
		nudge = ndsplineeval(table.get(), tablecoords, centers, 0);
		ENSURE(nudge >= 0);
		
		ENSURE_DISTANCE(base, nudge, std::numeric_limits<float>::epsilon(),
		    "Time quantile is continuous "
		    "(slope is less than FLOAT_EPSILON/DBL_EPSILON)");
		
		base = nudge;
		
		tablecoords[time_dim] = table->knots[time_dim][i]*(1+eps);
		ENSURE_EQUAL(tablesearchcenters(table.get(), tablecoords, centers), 0,
		    "tablesearchcenters() succeeds inside the knot field.");
		
		nudge = ndsplineeval(table.get(), tablecoords, centers, 0);
		ENSURE(nudge >= 0);
		
		ENSURE_DISTANCE(base, nudge, std::numeric_limits<float>::epsilon(),
		    "Time quantile is continuous "
		    "(slope is less than FLOAT_EPSILON/DBL_EPSILON)");
	}

	/* Check the transition into the knot field from the right */
	tablecoords[time_dim] = table->knots[time_dim][table->nknots[time_dim]-1]*(1+eps);
	ENSURE(tablecoords[time_dim] > table->knots[time_dim][table->nknots[time_dim]-1]);
	ENSURE(tablesearchcenters(table.get(), tablecoords, centers) != 0,
	    "tablesearchcenters() fails right at the right edge of the knot field.");
	
	tablecoords[time_dim] = table->knots[time_dim][table->nknots[time_dim]-1];
	
	ENSURE_EQUAL(tablesearchcenters(table.get(), tablecoords, centers), 0,
	    "tablesearchcenters() succeeds just inside the right edge of the knot field.");
	ENSURE_EQUAL(centers[time_dim],
	    table->nknots[time_dim]-table->order[time_dim]-2,
	    "centers[time_dim] holds the first fully-supported knot index.");
	
	base = ndsplineeval(table.get(), tablecoords, centers, 0);
	
	ENSURE(base >= 0, "Base quantile is positive.");
	ENSURE_DISTANCE(0, base, std::numeric_limits<float>::epsilon(),
	    "Time quantile is continuous "
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

/*
 * bspline_nonzero() can be made to return sensical values anywhere
 * in the knot field.
 */

TEST(bspline_nonzero_vs_bspline)
{
	unsigned i;
	const unsigned n_knots = 10;
	const int order = 2;
	double x, *knots;
	int center, offset;
	float localbasis_bsplvb[order+1], localbasis_bspline[order+1],
	    localbasis_bsplvb_deriv[order+1], localbasis_bspline_deriv[order+1];
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
	
		bspline_nonzero(knots, n_knots, x, center /* left */,
		    order /* jhigh */, localbasis_bsplvb, localbasis_bsplvb_deriv);
	
		for (offset = -order; offset <= 0; offset++) {
			ENSURE(offset+order >= 0);
			localbasis_bspline[offset+order] =
			    bspline(knots, x, center + offset, order);
			localbasis_bspline_deriv[offset+order] =
			    bspline_deriv(knots, x, center + offset, order);
		}
	
		for (offset = 0; offset < order+1; offset++) {
			ENSURE_DISTANCE(localbasis_bspline[offset], localbasis_bsplvb[offset], tol);
			ENSURE_DISTANCE(localbasis_bspline_deriv[offset], localbasis_bsplvb_deriv[offset], tol);
			
		}
	}
	
	/* Within the support. */
	for (i = order+1; i < n_knots-order-1; i++) {
		x = (knots[i]+knots[i+1])/2.0;
		center = i;
		
		bspline_nonzero(knots, n_knots, x, center /* left */,
		    order /* jhigh */, localbasis_bsplvb, localbasis_bsplvb_deriv);
	
		for (offset = -order; offset <= 0; offset++) {
			ENSURE(offset+order >= 0);
			localbasis_bspline[offset+order] =
			    bspline(knots, x, center + offset, order);
			localbasis_bspline_deriv[offset+order] =
			    bspline_deriv(knots, x, center + offset, order);
		}
	
		for (offset = 0; offset < order+1; offset++) {
			ENSURE_DISTANCE(localbasis_bspline[offset], localbasis_bsplvb[offset], tol);
			ENSURE_DISTANCE(localbasis_bspline_deriv[offset], localbasis_bsplvb_deriv[offset], tol);
			
		}
	}

	/* After the last first fully-supported knot. */
	for (i = n_knots-order-1; i < n_knots-2; i++) {
		x = (knots[i]+knots[i+1])/2.0; /* Less than fully-supported */
		center = n_knots-order-2; /* Last fully-supported spline. */
	
		ENSURE(int(i) >= center);
	
		bspline_nonzero(knots, n_knots, x, center /* left */,
		    order /* jhigh */, localbasis_bsplvb, localbasis_bsplvb_deriv);
	
		for (offset = -order; offset <= 0; offset++) {
			ENSURE(offset+order >= 0);
			localbasis_bspline[offset+order] =
			    bspline(knots, x, center + offset, order);
			localbasis_bspline_deriv[offset+order] =
			    bspline_deriv(knots, x, center + offset, order);
		}
	
		for (offset = 0; offset < order+1; offset++) {
			ENSURE_DISTANCE(localbasis_bspline[offset], localbasis_bsplvb[offset], tol);
			ENSURE_DISTANCE(localbasis_bspline_deriv[offset], localbasis_bsplvb_deriv[offset], tol);
			
		}
	}
}

