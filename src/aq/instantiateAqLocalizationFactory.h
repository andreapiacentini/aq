/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef AQ_MODEL_INSTANTIATEAQLOCALIZATIONFACTORY_H_
#define AQ_MODEL_INSTANTIATEAQLOCALIZATIONFACTORY_H_

#include "oops/interface/LocalizationBase.h"

#include "oops/aq/AqTraits.h"
#include "oops/aq/LocalizationMatrixAQ.h"

namespace aq {

void instantiateAqLocalizationFactory() {
  static oops::interface::LocalizationMaker<aq::AqTraits, LocalizationMatrixAQ> makerAQ_("AQ");
}

}  // namespace aq

#endif  // AQ_MODEL_INSTANTIATEAQLOCALIZATIONFACTORY_H_
