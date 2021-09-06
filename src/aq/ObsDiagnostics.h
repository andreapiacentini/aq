/*
 * (C) Copyright 2018  UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef AQ_OBSDIAGNOSTICS_H_
#define AQ_OBSDIAGNOSTICS_H_

#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/Printable.h"

#include "aq/ObsSpace.h"

namespace aq {
  class Locations;

// -----------------------------------------------------------------------------

class ObsDiagnostics : public util::Printable {
 public:
  ObsDiagnostics(const ObsSpace &, const Locations &, const oops::Variables &) {}
  ~ObsDiagnostics() {}

// I/O
  void save(const std::string &) const {}

 private:
  void print(std::ostream &) const {}
};
// -----------------------------------------------------------------------------
}  // namespace aq

#endif  // AQ_OBSDIAGNOSTICS_H_
