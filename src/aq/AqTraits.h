/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef AQ_AQTRAITS_H_
#define AQ_AQTRAITS_H_

#include <string>

#include "aq/AnalyticInit.h"
#include "aq/ErrorCovarianceAQ.h"
#include "aq/GeometryAQ.h"
#include "aq/GeometryAQIterator.h"
#include "aq/GetValuesAQ.h"
#include "aq/GetValuesTLAD.h"
#include "aq/GomAQ.h"
#include "aq/IncrementAQ.h"
#include "aq/LocationsAQ.h"
#include "aq/ModelBias.h"
#include "aq/ModelBiasCovariance.h"
#include "aq/ModelBiasIncrement.h"
#include "aq/ObsBias.h"
#include "aq/ObsBiasCovariance.h"
#include "aq/ObsBiasIncrement.h"
#include "aq/ObsDataAQ.h"
#include "aq/ObsDiagsAQ.h"
#include "aq/ObsIteratorAQ.h"
#include "aq/ObsOperatorAQ.h"
#include "aq/ObsOperatorTLAD.h"
#include "aq/ObsSpaceAQ.h"
#include "aq/ObsVecAQ.h"
#include "aq/StateAQ.h"

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

  typedef aq::ModelBias             ModelAuxControl;
  typedef aq::ModelBiasIncrement    ModelAuxIncrement;
  typedef aq::ModelBiasCovariance   ModelAuxCovariance;
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

#endif  // AQ_AQTRAITS_H_
