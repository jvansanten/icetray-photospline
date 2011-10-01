
#include <string>
extern "C" {
#include <photospline/splinetable.h>
}

class I3SplineTable {
public:
	I3SplineTable(const std::string &path);
	virtual ~I3SplineTable() {};
	int Eval(double *x, double *result);
protected:
	boost::shared_ptr<splinetable> table_;	
};
