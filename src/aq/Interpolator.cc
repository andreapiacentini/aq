/*
 * (C) Copyright 2022-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "aq/Interpolator.h"

#include <ostream>
#include <vector>

#include "aq/aq_fields_interface.h"
#include "aq/Geometry.h"
#include "aq/Increment.h"
#include "aq/State.h"

namespace aq {

// -----------------------------------------------------------------------------

Interpolator::Interpolator(const eckit::Configuration &, const Geometry & grid,
                           const std::vector<double> & lats, const std::vector<double> & lons)
  : nlevs_(1), nlocs_(lats.size()), locs_(2 * nlocs_)
{
  ASSERT(lats.size() == lons.size());
  for (size_t jj = 0; jj < nlocs_; ++jj) {
    locs_[2 * jj] = lats[jj];
    locs_[2 * jj + 1] = lons[jj];
  }
}

// -----------------------------------------------------------------------------

Interpolator::~Interpolator() {}

// -----------------------------------------------------------------------------

void Interpolator::apply(const oops::Variables & vars, const State & xx,
                         const std::vector<bool> & mask,
                         std::vector<double> & values) const {
  const size_t nvals = vars.size() * nlevs_ * nlocs_;
  values.resize(nvals);
  aq_fields_getvals_f90(xx.fields().toFortran(), vars, nlocs_, locs_[0], nvals, values[0]);
}

// -----------------------------------------------------------------------------

void Interpolator::apply(const oops::Variables & vars, const Increment & dx,
                         const std::vector<bool> & mask,
                         std::vector<double> & values) const {
  const size_t nvals = vars.size() * nlevs_ * nlocs_;
  values.resize(nvals);
  aq_fields_getvals_f90(dx.fields().toFortran(), vars, nlocs_, locs_[0], nvals, values[0]);
}

// -----------------------------------------------------------------------------

void Interpolator::applyAD(const oops::Variables & vars, Increment & dx,
                           const std::vector<bool> & mask,
                           const std::vector<double> & values) const {
  const size_t nvals = vars.size() * nlevs_ * nlocs_;
  ASSERT(values.size() == nvals);
  aq_fields_getvalsad_f90(dx.fields().toFortran(), vars, nlocs_, locs_[0], nvals, values[0]);
}

// -----------------------------------------------------------------------------

void Interpolator::print(std::ostream & os) const {
  os << "Interpolator";
}

// -----------------------------------------------------------------------------

}  // namespace aq
