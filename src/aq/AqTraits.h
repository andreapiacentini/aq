/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef AQ_MODEL_AQTRAITS_H_
#define AQ_MODEL_AQTRAITS_H_

#include <string>

#include "oops/aq/AnalyticInit.h"
#include "oops/aq/ErrorCovarianceAQ.h"
#include "oops/aq/GeometryAQ.h"
#include "oops/aq/GeometryAQIterator.h"
#include "oops/aq/GetValuesAQ.h"
#include "oops/aq/GetValuesTLAD.h"
#include "oops/aq/GomAQ.h"
#include "oops/aq/IncrementAQ.h"
#include "oops/aq/LocationsAQ.h"
#include "oops/aq/ObsBias.h"
#include "oops/aq/ObsBiasCovariance.h"
#include "oops/aq/ObsBiasIncrement.h"
#include "oops/aq/ObsDataAQ.h"
#include "oops/aq/ObsDiagsAQ.h"
#include "oops/aq/ObsIteratorAQ.h"
#include "oops/aq/ObsOperatorAQ.h"
#include "oops/aq/ObsOperatorTLAD.h"
#include "oops/aq/ObsSpaceAQ.h"
#include "oops/aq/ObsVecAQ.h"
#include "oops/aq/StateAQ.h"

namespace aq {

struct AqTraits {
  static std::string name() {return "AQ";}
  static std::string nameCovar() {return "AqError";}
  static std::string nameCovar4D() {return "AqError";}

  typedef aq::GeometryAQ            Geometry;

  typedef aq::GeometryAQIterator    GeometryIterator;

  typedef aq::GetValuesAQ           GetValues;
  typedef aq::GetValuesTLAD         LinearGetValues;

  typedef aq::StateAQ               State;
  typedef aq::IncrementAQ           Increment;
  typedef aq::ErrorCovarianceAQ     Covariance;
};

struct AqObsTraits {
  static std::string name() {return "AQ obs";}

  typedef aq::ObsSpaceAQ            ObsSpace;
  typedef aq::ObsVecAQ              ObsVector;
  typedef aq::ObsOperatorAQ         ObsOperator;
  typedef aq::ObsOperatorTLAD       LinearObsOperator;
  template <typename DATATYPE> using ObsDataVector = aq::ObsDataAQ<DATATYPE>;
  typedef aq::ObsIteratorAQ         GeometryIterator;

  typedef aq::ObsBias               ObsAuxControl;
  typedef aq::ObsBiasIncrement      ObsAuxIncrement;
  typedef aq::ObsBiasCovariance     ObsAuxCovariance;

  typedef aq::ObsDiagsAQ            ObsDiagnostics;

  typedef aq::GomAQ                 GeoVaLs;
  typedef aq::LocationsAQ           Locations;
  typedef aq::AnalyticInit          AnalyticInit;
};

}  // namespace aq

#endif  // AQ_MODEL_AQTRAITS_H_
