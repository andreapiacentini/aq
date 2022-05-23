/*
 * This file is part of the Air Quality Ensemble Data Assimilation suite AQ.
 *
 * (C) Copyright 2022 CERFACS
 *
 * AQ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.
 *
 * AQ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * A copy of the GNU Lesser General Public License is distributed
 * along with AQ (files LICENSE.md, COPYING and COPYING.LESSER).
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
  typedef aq::Locations            Locations;
};

}  // namespace aq

#endif  // AQ_TRAITSFWD_H_
