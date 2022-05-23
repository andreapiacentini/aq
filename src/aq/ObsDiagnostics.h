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
