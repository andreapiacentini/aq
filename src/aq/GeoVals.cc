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

#include "aq/GeoVals.h"

#include <iomanip>

#include "aq/aq_geovals_interface.h"
#include "aq/Locations.h"
#include "aq/ObsSpace.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

namespace aq {

// -----------------------------------------------------------------------------
GeoVals::GeoVals(const Locations & locs, const oops::Variables & vars,
             const std::vector<size_t> & sizes):
  vars_(vars), comm_(locs.comm())
{
  // geovals_setup just creates and allocates the GeoVaLs object without filling
  // in values
  aq_geovals_setup_f90(keyGeoVals_, locs, vars_, &comm_);
}
// -----------------------------------------------------------------------------
/*! AQ GeoVaLs Constructor with Config */

  GeoVals::GeoVals(const Parameters_ & params,
               const ObsSpace & ospace, const oops::Variables & vars):
  vars_(vars), comm_(ospace.comm())
{
  aq_geovals_create_f90(keyGeoVals_, vars_, &comm_);
  aq_geovals_read_file_f90(keyGeoVals_, params.toConfiguration());
}
// -----------------------------------------------------------------------------
// Copy constructor
GeoVals::GeoVals(const GeoVals & other):
  vars_(other.vars_), comm_(other.comm_)
{
  aq_geovals_create_f90(keyGeoVals_, vars_, &comm_);
  aq_geovals_copy_f90(keyGeoVals_, other.keyGeoVals_);
}
// -----------------------------------------------------------------------------
GeoVals::~GeoVals() {
  aq_geovals_delete_f90(keyGeoVals_);
}
// -----------------------------------------------------------------------------
double GeoVals::rms() const {
  double zz;
  aq_geovals_rms_f90(keyGeoVals_, zz);
  return zz;
}
// -----------------------------------------------------------------------------
double GeoVals::normalizedrms(const GeoVals & rhs) const {
  double zz;
  GeoVals temp_GeoVals(*this);
  aq_geovals_divide_f90(temp_GeoVals.keyGeoVals_, rhs.keyGeoVals_);
  zz = temp_GeoVals.rms();
  return zz;
}
// -----------------------------------------------------------------------------
void GeoVals::zero() {
  aq_geovals_zero_f90(keyGeoVals_);
}
// -----------------------------------------------------------------------------
void GeoVals::random() {
  aq_geovals_random_f90(keyGeoVals_);
}
// -----------------------------------------------------------------------------
GeoVals & GeoVals::operator=(const GeoVals & rhs) {
  const int keyGeoValsRhs = rhs.keyGeoVals_;
  aq_geovals_copy_f90(keyGeoVals_, keyGeoValsRhs);
  return *this;
}
// -----------------------------------------------------------------------------
GeoVals & GeoVals::operator*=(const double & zz) {
  aq_geovals_mult_f90(keyGeoVals_, zz);
  return *this;
}
// -----------------------------------------------------------------------------
GeoVals & GeoVals::operator+=(const GeoVals & other) {
  aq_geovals_add_f90(keyGeoVals_, other.keyGeoVals_);
  return *this;
}
// -----------------------------------------------------------------------------
GeoVals & GeoVals::operator-=(const GeoVals & other) {
  aq_geovals_diff_f90(keyGeoVals_, other.keyGeoVals_);
  return *this;
}
// -----------------------------------------------------------------------------
GeoVals & GeoVals::operator*=(const GeoVals & other) {
  aq_geovals_schurmult_f90(keyGeoVals_, other.keyGeoVals_);
  return *this;
}
// -----------------------------------------------------------------------------
double GeoVals::dot_product_with(const GeoVals & other) const {
  double zz;
  aq_geovals_dotprod_f90(keyGeoVals_, other.keyGeoVals_, zz);
  return zz;
}

// -----------------------------------------------------------------------------
void GeoVals::fill(const std::vector<size_t> & indx, const std::vector<double> & vals) {
  const size_t npts = indx.size();
  const size_t nvals = vals.size();
  std::vector<int> findx(indx.size());
  for (size_t jj = 0; jj < indx.size(); ++jj) findx[jj] = indx[jj] + 1;
  aq_geovals_fill_f90(keyGeoVals_, npts, findx[0], nvals, vals[0]);
}
// -----------------------------------------------------------------------------
void GeoVals::fillAD(const std::vector<size_t> & indx, std::vector<double> & vals) const {
  const size_t npts = indx.size();
  const size_t nvals = vals.size();
  std::vector<int> findx(indx.size());
  for (size_t jj = 0; jj < indx.size(); ++jj) findx[jj] = indx[jj] + 1;
  aq_geovals_fillad_f90(keyGeoVals_, npts, findx[0], nvals, vals[0]);
}
// -----------------------------------------------------------------------------
void GeoVals::read(const Parameters_ & params) {
  if (comm_.rank() == 0) {
    aq_geovals_read_file_f90(keyGeoVals_, params.toConfiguration());
  }
}
// -----------------------------------------------------------------------------
void GeoVals::write(const Parameters_ & params) const {
  if (comm_.rank() == 0) {
    aq_geovals_write_file_f90(keyGeoVals_, params.toConfiguration());
  }
}
// -----------------------------------------------------------------------------
void GeoVals::print(std::ostream & os) const {
  int nobs;
  double zmin, zmax, zave, zstd;
  aq_geovals_stats_f90(keyGeoVals_, nobs, zmin, zmax, zave, zstd);
  std::ios_base::fmtflags f(os.flags());
  os << "  GeoVals for var " << vars_[0] << ": nobs= " << nobs
     << std::scientific << std::setprecision(4)
     << "  Min=" << std::setw(12) << zmin
     << ", Max=" << std::setw(12) << zmax
     << ", Ave=" << std::setw(12) << zave
     << ", Std=" << std::setw(12) << zstd;
  os.flags(f);

  // If the min value across all variables is positive, then this may be an
  // error measurement.  If so, print the location and variable where the
  // maximum occurs to the debug stream, for use in debugging

  if (zmin >= 0.0) {
    double mxval;
    int iloc;
    oops::Variables maxvar;

    aq_geovals_maxloc_f90(keyGeoVals_, mxval, iloc, maxvar);

    oops::Log::debug() << "GeoVals: Maximum Value = " << std::setprecision(4)
                       << mxval << " at location = " << iloc
                       << " and variable = " << maxvar << std::endl;
  }
}
// -----------------------------------------------------------------------------
}  // namespace aq
