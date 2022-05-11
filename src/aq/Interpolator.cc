/*
 * (C) Copyright 2022-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "aq/Interpolator.h"

#include <ostream>
#include <vector>

#include "aq/aq_interpolator_interface.h"
#include "aq/Geometry.h"
#include "aq/Increment.h"
#include "aq/State.h"

namespace aq {

// -----------------------------------------------------------------------------

Interpolator::Interpolator(const eckit::Configuration & conf, const Geometry & grid,
                           const std::vector<double> & lats, const std::vector<double> & lons)
  : nlevs_(1), geom_(new Geometry(grid)), nlocs_(lats.size()), lats_(lats), lons_(lons)
{
  ASSERT(lats.size() == lons.size());
  aq_interpolator_create_f90(keyInterp_, conf, geom_->toFortran(), nlocs_, lats_[0], lons_[0]);
}

// -----------------------------------------------------------------------------

Interpolator::~Interpolator() {
  aq_interpolator_delete_f90(keyInterp_);
}

// -----------------------------------------------------------------------------

void Interpolator::apply(const oops::Variables & vars, const State & xx,
                         const std::vector<bool> & mask,
                         std::vector<double> & values) const {
  const size_t nvals = vars.size() * nlevs_ * nlocs_;
  values.resize(nvals);
  ASSERT(mask.size() == values.size());
  std::vector<int> imask(nvals, 0);
  for (size_t jobs = 0; jobs < nvals; ++jobs) {
    if (mask[jobs]) imask[jobs] = 1;
  }
  aq_interpolator_apply_f90(keyInterp_, xx.fields().toFortran(), vars,
                            nvals, imask[0], values[0]);
}

// -----------------------------------------------------------------------------

void Interpolator::apply(const oops::Variables & vars, const Increment & dx,
                         const std::vector<bool> & mask,
                         std::vector<double> & values) const {
  const size_t nvals = vars.size() * nlevs_ * nlocs_;
  values.resize(nvals);
  ASSERT(mask.size() == values.size());
  std::vector<int> imask(nvals, 0);
  for (size_t jobs = 0; jobs < nvals; ++jobs) {
    if (mask[jobs]) imask[jobs] = 1;
  }
  aq_interpolator_apply_f90(keyInterp_, dx.fields().toFortran(), vars,
                            nvals, imask[0], values[0]);
}

// -----------------------------------------------------------------------------

void Interpolator::applyAD(const oops::Variables & vars, Increment & dx,
                           const std::vector<bool> & mask,
                           const std::vector<double> & values) const {
  const size_t nvals = vars.size() * nlevs_ * nlocs_;
  ASSERT(values.size() == nvals);
  ASSERT(mask.size() == values.size());
  std::vector<int> imask(nvals, 0);
  for (size_t jobs = 0; jobs < nvals; ++jobs) {
    if (mask[jobs]) imask[jobs] = 1;
  }
  aq_interpolator_applyAD_f90(keyInterp_, dx.fields().toFortran(), vars,
                              nvals, imask[0], values[0]);
}

// -----------------------------------------------------------------------------

void Interpolator::print(std::ostream & os) const {
  os << "Interpolator";
}

// -----------------------------------------------------------------------------

}  // namespace aq
