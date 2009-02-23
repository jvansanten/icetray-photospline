/**
 *@brief Wrappers around PSInterface classes to define a python module to acces Photoncis tables
 *
 * copyright  (C) 2007
 * the icecube collaboration
 *
 * $Id: PyPSI.cxx 45789 2008-05-30 10:14:42Z bvoigt $
 *
 * @version $Revision: 45789 $
 * @date $LastChangedDate: 2008-05-30 12:14:42 +0200 (Fri, 30 May 2008) $
 * @author Bernhard Voigt <bernhard.voigt@desy.de>   Last changed by: $LastChangedBy: bvoigt $
 */

namespace bp=boost::python;

#include "PSInterface/PSI_Coordinate.h"
#include "PSInterface/PSI_Coordinate_Photonics.h"
#include "PSInterface/PSInterface.h"
#include "PSInterface/PSI_Photonics.h"
#include "PSInterface/PSI_PhotonicsPhotorec.h"


/* I couldn't figure out how to deal with 'pass by reference values' of methods.
 * Hence I came up with simple wrapper functions that take as a first argument
 * an instance of the actual object that defines the method to be wrapped. Other
 * arguments being no references are just passed through.
 *
 * Since object method's in python are called with the first argument being an object
 * instance, this works without problems.
 */

// PSI_Coordinate method wrappers
bp::tuple CalcTrackDirection(PSI_Coordinate& coordinate) {
  double x, y, z = 0;
  coordinate.CalcTrackDirection(x, y, z);
  return bp::make_tuple(x, y, z);
}

bp::tuple CalcRhoVector(PSI_Coordinate& coordinate) {
  double x, y, z = 0;
  coordinate.CalcRhoVector(x, y, z);
  return bp::make_tuple(x, y, z);
}

bp::tuple CalculatePPerp(PSI_Coordinate_Photonics& coordinate) {
  double x, y, z = 0;
  coordinate.CalculatePPerp(x, y, z);
  return bp::make_tuple(x, y, z);
}

// PSI_Photonics method wrappers
PSI_Coordinate* MakeCoordinate(PSI_Photonics& psi, const double& x, const double& y, const double& z, const double& o,
			       const double& tx, const double& ty, const double& tz, const double& tze,
			       const double& taz, const double& tl, const double& te, const int& tt) {
  return psi.MakeCoordinate(x, y, z, o, tx, ty, tz, tze, taz, tl, te, tt);
}

double GetMeanAmplitude(PSI_Photonics& psi, PSI_Coordinate* coordinate) {
  double amp = -1.0;
  psi.GetMeanAmplitude(amp, coordinate);
  return amp;
}

double GetTimeDelay(PSI_Photonics& psi, const double& random, PSI_Coordinate* coordinate) {
  double delay = -1.0;
  psi.GetTimeDelay(delay, random, coordinate);
  return delay;
}

double GetProbability(PSI_Photonics& psi, const double& delay, PSI_Coordinate* coordinate) {
  double prob = -1.0;
  psi.GetProbability(prob, delay, coordinate);
  return prob;
}

bp::tuple GetRefRefractive(PSI_Photonics& psi, const int& level) {
  double ref_ng, ref_np = 0;
  psi.GetRefRefractive(level, ref_ng, ref_np);
  // return a tuple of the two values asked for
  return bp::make_tuple(ref_ng, ref_np);
}

bp::tuple GetAngularInterval(PSI_Photonics& psi, const int& level) {
  float low, high = 0;
  psi.GetAngularInterval(level, low, high);
  // return a tuple of the two values asked for
  return bp::make_tuple(low, high);
}

bp::tuple GetDepthInterval(PSI_Photonics& psi, const int& level) {
  float low, high = 0;
  psi.GetDepthInterval(level, low, high);
  // return a tuple of the two values asked for
  return bp::make_tuple(low, high);
}

// Photorec
bp::tuple Get(PSI_PhotonicsPhotorec& photorec, PSI_Coordinate_Photonics* coordinate, const double& delay) {
    double amplitude = -1;
    double probability = -1;

    photorec.Get(coordinate, delay, amplitude, probability);

    return bp::make_tuple(amplitude, probability);
}


/* Define a python module named libpycmc
 * 
 * The following classes with public member functions are defined:
 * PSI_Coordinate
 * PSI_Coordinate_Photonics
 * PSI_Photonics
 */
BOOST_PYTHON_MODULE(PSInterface) {

  bp::class_<PSI_Coordinate >("PSI_Coordinate", 
			      bp::init<double, double, double, double, 
			      double, double, double, double, double, double, double, int >())
    .def("Clear", &PSI_Coordinate::Clear)
    .def("Set", &PSI_Coordinate::Set)
    .def("CalcTrackDirection", &CalcTrackDirection)  // use the wrapper function to deal with pass by reference args
    .def("CalcDistance", &PSI_Coordinate::CalcDistance)
    .def("CalcOnTrackDistance", &PSI_Coordinate::CalcOnTrackDistance)
    .def("CalcRho", &PSI_Coordinate::CalcRho)
    .def("CalcRhoVector", &CalcRhoVector)  // use the wrapper function to deal with pass by reference args
    .def("CalcTrackToOpticalModuleVector", &PSI_Coordinate::CalcTrackToOpticalModuleVector)
    .def("CalcTGeo", &PSI_Coordinate::CalcTGeo)
    .def("GetOpticalModuleX", &PSI_Coordinate::GetOpticalModuleX)
    .def("GetOpticalModuleY", &PSI_Coordinate::GetOpticalModuleY)
    .def("GetOpticalModuleZ", &PSI_Coordinate::GetOpticalModuleZ)
    .def("GetOpticalModuleOrientation", &PSI_Coordinate::GetOpticalModuleOrientation)
    .def("GetTrackX", &PSI_Coordinate::GetTrackX)
    .def("GetTrackY", &PSI_Coordinate::GetTrackY)
    .def("GetTrackZ", &PSI_Coordinate::GetTrackZ)
    .def("GetTrackTheta", &PSI_Coordinate::GetTrackTheta)
    .def("GetTrackPhi", &PSI_Coordinate::GetTrackPhi)
    .def("GetTrackLength", &PSI_Coordinate::GetTrackLength)
    .def("GetTrackEnergy", &PSI_Coordinate::GetTrackEnergy)
    .def("GetTrackType", &PSI_Coordinate::GetTrackType)
    .def("MuonEnergyCorrection", &PSI_Coordinate::MuonEnergyCorrection);

  // derive class of PSI_Coordinate, use boost.python bases class for that
  bp::class_<PSI_Coordinate_Photonics, bp::bases<PSI_Coordinate> >("PSI_Coordinate_Photonics",   
								   bp::init<double, double, double, double, 
								   double, double, double, double, double, double, double, int >())
    
    .def("CalculatePPerp", &CalculatePPerp)  // use the wrapper function to deal with pass by reference args
    .def("CalculateCoordinate", &PSI_Coordinate_Photonics::CalculateCoordinate)
    .def("GetCoordinateSrc", &PSI_Coordinate_Photonics::GetCoordinateZSrc)
    .def("GetCoordinateTheta", &PSI_Coordinate_Photonics::GetCoordinateTheta)
    .def("GetCoordinateL", &PSI_Coordinate_Photonics::GetCoordinateL)
    .def("GetCoordinateRho", &PSI_Coordinate_Photonics::GetCoordinateRho)
    .def("GetCoordinatePhi", &PSI_Coordinate_Photonics::GetCoordinatePhi)
    .def("GetCoordinateDist", &PSI_Coordinate_Photonics::GetCoordinateDist)
    .def("GetCoordinateStopZSrc", &PSI_Coordinate_Photonics::GetCoordinateStopZSrc)
    .def("GetCoordinateStopL", &PSI_Coordinate_Photonics::GetCoordinateStopL)
    .def("GetCoordinateStopDist", &PSI_Coordinate_Photonics::GetCoordinateStopDist)
    .def("GetZOrigin", &PSI_Coordinate_Photonics::GetZOrigin)
    .def("IsFinite", &PSI_Coordinate_Photonics::IsFinite)
    .def("SetZOrigin", &PSI_Coordinate_Photonics::SetZOrigin);
  

  // boost::noncopyable tells python that these instances can't be copied, as it is the case since 
  // the copy constructor is private and this gives a compiler error anyway
  bp::class_<PSI_Photonics, boost::noncopyable>("PSI_Photonics")
    .def("MakeCoordinate", &MakeCoordinate, bp::return_value_policy<bp::manage_new_object>())
    .def("GetMeanAmplitude", &GetMeanAmplitude)  // use the wrapper function to deal with pass by reference args
    .def("GetTimeDelay", &GetTimeDelay)		 // use the wrapper function to deal with pass by reference args
    .def("GetProbability", &GetProbability)	 // use the wrapper function to deal with pass by reference args
    .def("LoadTables", &PSI_Photonics::LoadTables)
    .def("ClearTables", &PSI_Photonics::ClearTables)
    .def("GetMemoryUsage", &PSI_Photonics::GetMemoryUsage)
    .def("SetAngularSelection", &PSI_Photonics::SetAngularSelection)
    .def("SetDepthSelection", &PSI_Photonics::SetDepthSelection)
    .def("GetAngularInterval", &GetAngularInterval)  // use the wrapper function to deal with pass by reference args
    .def("GetDepthInterval", &GetDepthInterval)	     // use the wrapper function to deal with pass by reference args
    .def("GetAngularMemoryRequirement", &PSI_Photonics::GetAngularMemoryRequirement)
    .def("SetPhotonicsVerbosity", &PSI_Photonics::SetPhotonicsVerbosity)
    .def("SetRefRefractive", &PSI_Photonics::SetRefRefractive)
    .def("GetRefRefractive", &GetRefRefractive)      // use the wrapper function to deal with pass by reference args
    .def("IsTablesLoaded", &PSI_Photonics::IsTablesLoaded)
    .def("SetDriverFileName", &PSI_Photonics::SetDriverFileName)
    .def("GetDriverFileName", &PSI_Photonics::GetDriverFileName)
    .def("GetInterpolationMode", &PSI_Photonics::GetInterpolationMode)
    .def("SetInterpolationMode", &PSI_Photonics::SetInterpolationMode);

  bp::class_<PSI_PhotonicsPhotorec, bp::bases<PSI_Photonics>, boost::noncopyable>("PSI_PhotonicsPhotorec")
    .def("Get", &Get);
}
