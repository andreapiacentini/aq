/*
 * (C) Copyright 2018  UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef AQ_MODEL_OBSDIAGSAQ_H_
#define AQ_MODEL_OBSDIAGSAQ_H_

#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/Printable.h"

#include "oops/aq/ObsSpaceAQ.h"

namespace aq {
  class LocationsAQ;

// -----------------------------------------------------------------------------

class ObsDiagsAQ : public util::Printable {
 public:
  ObsDiagsAQ(const ObsSpaceAQ &, const LocationsAQ &, const oops::Variables &) {}
  ~ObsDiagsAQ() {}

// I/O
  void save(const std::string &) const {}

 private:
  void print(std::ostream &) const {}
};
// -----------------------------------------------------------------------------
}  // namespace aq

#endif  // AQ_MODEL_OBSDIAGSAQ_H_
