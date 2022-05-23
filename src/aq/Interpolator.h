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

#ifndef AQ_INTERPOLATOR_H_
#define AQ_INTERPOLATOR_H_

#include <memory>
#include <ostream>
#include <vector>

#include "oops/util/Printable.h"

#include "aq/Geometry.h"
#include "aq/interface.h"

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
// Utilities
  const int & toFortran() const {return keyInterp_;}

 private:
  void print(std::ostream &) const;
  F90interp keyInterp_;
  std::shared_ptr<const Geometry> geom_;
  const size_t nlevs_;
  const size_t nlocs_;
  std::vector<double> lats_;
  std::vector<double> lons_;
};
// -----------------------------------------------------------------------------

}  // namespace aq
#endif  // AQ_INTERPOLATOR_H_
