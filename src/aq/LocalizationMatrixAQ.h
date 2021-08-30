/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2020-2020 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef AQ_MODEL_LOCALIZATIONMATRIXAQ_H_
#define AQ_MODEL_LOCALIZATIONMATRIXAQ_H_

#include <ostream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"

#include "oops/interface/LocalizationBase.h"

#include "oops/aq/AqFortran.h"
#include "oops/aq/AqTraits.h"
#include "oops/aq/GeometryAQ.h"

// Forward declarations
namespace aq {
  class GeometryAQ;
  class IncrementAQ;

/// Localization matrix for AQ model.

// -----------------------------------------------------------------------------
class LocalizationMatrixAQ: public oops::interface::LocalizationBase<aq::AqTraits> {
 public:
  static const std::string classname() {return "aq::LocalizationMatrixAQ";}

  LocalizationMatrixAQ(const GeometryAQ &, const eckit::Configuration &);
  ~LocalizationMatrixAQ();

  void randomize(IncrementAQ &) const override;
  void multiply(IncrementAQ &) const override;

 private:
  void print(std::ostream &) const override;
  F90lclz keyLocal_;
};
// -----------------------------------------------------------------------------

}  // namespace aq

#endif  // AQ_MODEL_LOCALIZATIONMATRIXAQ_H_
