/*
 * (C) Copyright 2022-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <ostream>
#include <vector>

#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace aq {
  class Geometry;
  class Increment;
  class State;

// -----------------------------------------------------------------------------

class Interpolator : public util::Printable {
 public:
  Interpolator(const eckit::Configuration &, const Geometry &,
               const std::vector<double> &, const std::vector<double> &);
  ~Interpolator();

  void apply(const oops::Variables &, const State &,
             const std::vector<bool> &, std::vector<double> &) const;
  void apply(const oops::Variables &, const Increment &,
             const std::vector<bool> &, std::vector<double> &) const;
  void applyAD(const oops::Variables &, Increment &,
               const std::vector<bool> &, const std::vector<double> &) const;

 private:
  void print(std::ostream &) const;

  const size_t nlevs_;
  const size_t nlocs_;
  std::vector<double> locs_;
};
// -----------------------------------------------------------------------------

}  // namespace aq
