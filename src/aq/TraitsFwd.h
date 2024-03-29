/*
 * (C) Copyright 2021 UCAR.
 * (C) Copyright 2021-2022 CERFACS.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef AQ_TRAITSFWD_H_
#define AQ_TRAITSFWD_H_

#include <string>

namespace aq {

class Geometry;
class GeometryIterator;

class State;
class Increment;
class ErrorCovariance;
class InterpolatorQG;

class ChangeVar;
class ChangeVarTLAD;

class ModelBias;
class ModelBiasIncrement;
class ModelBiasCovariance;

class ObsSpace;
class ObsVec;
template <typename DATATYPE> class ObsData;
class ObsIterator;

class ObsOperator;
class ObsOperatorTLAD;
class ObsBias;
class ObsBiasIncrement;
class ObsBiasCovariance;
class ObsBiasPreconditioner;
class ObsDiagnostics;

class GeoVals;
class Locations;

struct Traits {
  static std::string name() {return "AQ";}
  static std::string nameCovar() {return "AqError";}
  static std::string nameCovar4D() {return "AqError";}

  typedef aq::Geometry            Geometry;

  typedef aq::GeometryIterator    GeometryIterator;

  typedef aq::ChangeVar           VariableChange;
  typedef aq::LinearChangeVar     LinearVariableChange;

  typedef aq::State               State;
  typedef aq::Increment           Increment;
  typedef aq::Covariance          Covariance;
  typedef aq::Interpolator        LocalInterpolator;

  typedef aq::ModelAuxControl     ModelAuxControl;
  typedef aq::ModelAuxIncrement   ModelAuxIncrement;
  typedef aq::ModelAuxCovariance  ModelAuxCovariance;
  typedef aq::ModelData           ModelData;
};

struct ObsTraits {
  static std::string name() {return "AQ obs";}

  typedef aq::ObsSpace             ObsSpace;
  typedef aq::ObsVec               ObsVector;
  typedef aq::ObsOperator          ObsOperator;
  typedef aq::LinearObsOperator    LinearObsOperator;
  template <typename DATATYPE> using ObsDataVector = aq::ObsData<DATATYPE>;
  typedef aq::ObsIterator          GeometryIterator;

  typedef aq::ObsAuxControl        ObsAuxControl;
  typedef aq::ObsAuxIncrement      ObsAuxIncrement;
  typedef aq::ObsAuxCovariance     ObsAuxCovariance;
  typedef aq::ObsAuxPreconditioner ObsAuxPreconditioner;

  typedef aq::ObsDiagnostics       ObsDiagnostics;

  typedef aq::GeoVals              GeoVaLs;
  typedef aq::Locations            SampledLocations;
};

}  // namespace aq

#endif  // AQ_TRAITSFWD_H_
