/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef AQ_TRAITS_H_
#define AQ_TRAITS_H_

#include <string>

#include "aq/Covariance.h"
#include "aq/Geometry.h"
#include "aq/GeometryIterator.h"
#include "aq/GeoVals.h"
#include "aq/GetValues.h"
#include "aq/Increment.h"
#include "aq/LinearGetValues.h"
#include "aq/LinearObsOperator.h"
#include "aq/Locations.h"
#include "aq/ModelAuxControl.h"
#include "aq/ModelAuxCovariance.h"
#include "aq/ModelAuxIncrement.h"
#include "aq/ObsAuxControl.h"
#include "aq/ObsAuxCovariance.h"
#include "aq/ObsAuxIncrement.h"
#include "aq/ObsData.h"
#include "aq/ObsDiagnostics.h"
#include "aq/ObsIterator.h"
#include "aq/ObsOperator.h"
#include "aq/ObsSpace.h"
#include "aq/ObsVec.h"
#include "aq/State.h"

namespace aq {

struct Traits {
  static std::string name() {return "AQ";}
  static std::string nameCovar() {return "AqError";}
  static std::string nameCovar4D() {return "AqError";}

  typedef aq::Geometry            Geometry;

  typedef aq::GeometryIterator    GeometryIterator;

  typedef aq::GetValues           GetValues;
  typedef aq::LinearGetValues     LinearGetValues;

  typedef aq::State               State;
  typedef aq::Increment           Increment;
  typedef aq::Covariance     Covariance;

  typedef aq::ModelAuxControl             ModelAuxControl;
  typedef aq::ModelAuxIncrement    ModelAuxIncrement;
  typedef aq::ModelAuxCovariance   ModelAuxCovariance;
};

struct ObsTraits {
  static std::string name() {return "AQ obs";}

  typedef aq::ObsSpace            ObsSpace;
  typedef aq::ObsVec              ObsVector;
  typedef aq::ObsOperator         ObsOperator;
  typedef aq::LinearObsOperator       LinearObsOperator;
  template <typename DATATYPE> using ObsDataVector = aq::ObsData<DATATYPE>;
  typedef aq::ObsIterator         GeometryIterator;

  typedef aq::ObsAuxControl               ObsAuxControl;
  typedef aq::ObsAuxIncrement      ObsAuxIncrement;
  typedef aq::ObsAuxCovariance     ObsAuxCovariance;

  typedef aq::ObsDiagnostics            ObsDiagnostics;

  typedef aq::GeoVals               GeoVaLs;
  typedef aq::Locations           Locations;
};

}  // namespace aq

#endif  // AQ_TRAITS_H_
