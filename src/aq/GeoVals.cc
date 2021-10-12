/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "aq/GeoVals.h"

#include <iomanip>

#include "aq/aq_geovals_interface.h"
#include "aq/Locations.h"
#include "aq/ObsSpace.h"
#include "eckit/config/Configuration.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

namespace aq {

// -----------------------------------------------------------------------------
GeoVals::GeoVals(const Locations & locs, const oops::Variables & vars,
             const std::vector<size_t> & sizes):
  vars_(vars), comm_(locs.comm())
{
  if (comm_.rank() == 0) {
    // geovals_setup just creates and allocates the GeoVaLs object without filling
    // in values
    aq_geovals_setup_f90(keyGeoVals_, locs, vars_);
  }
}
// -----------------------------------------------------------------------------
/*! AQ GeoVaLs Constructor with Config */

  GeoVals::GeoVals(const eckit::Configuration & config,
               const ObsSpace & ospace, const oops::Variables & vars):
  vars_(vars), comm_(ospace.comm())
{
  if (comm_.rank() == 0) {
    aq_geovals_create_f90(keyGeoVals_, vars_);
    aq_geovals_read_file_f90(keyGeoVals_, config);
  }
}
// -----------------------------------------------------------------------------
// Copy constructor
GeoVals::GeoVals(const GeoVals & other):
  vars_(other.vars_), comm_(other.comm_)
{
  if (comm_.rank() == 0) {
    aq_geovals_create_f90(keyGeoVals_, vars_);
    aq_geovals_copy_f90(keyGeoVals_, other.keyGeoVals_);
  }
}
// -----------------------------------------------------------------------------
GeoVals::~GeoVals() {
  if (comm_.rank() == 0) {
    aq_geovals_delete_f90(keyGeoVals_);
  }
}
// -----------------------------------------------------------------------------
double GeoVals::rms() const {
  double zz;
  if (comm_.rank() == 0) {
    aq_geovals_rms_f90(keyGeoVals_, zz);
  }
  comm_.broadcast(zz, 0);
  return zz;
}
// -----------------------------------------------------------------------------
double GeoVals::normalizedrms(const GeoVals & rhs) const {
  double zz;
  GeoVals temp_GeoVals(*this);
  if (comm_.rank() == 0) {
    aq_geovals_divide_f90(temp_GeoVals.keyGeoVals_, rhs.keyGeoVals_);
  }
  zz = temp_GeoVals.rms();
  return zz;
}
// -----------------------------------------------------------------------------
void GeoVals::zero() {
  if (comm_.rank() == 0) {
    aq_geovals_zero_f90(keyGeoVals_);
  }
}
// -----------------------------------------------------------------------------
void GeoVals::random() {
  if (comm_.rank() == 0) {
    aq_geovals_random_f90(keyGeoVals_);
  }
}
// -----------------------------------------------------------------------------
GeoVals & GeoVals::operator=(const GeoVals & rhs) {
  if (comm_.rank() == 0) {
    const int keyGeoValsRhs = rhs.keyGeoVals_;
    aq_geovals_copy_f90(keyGeoVals_, keyGeoValsRhs);
  }
  return *this;
}
// -----------------------------------------------------------------------------
GeoVals & GeoVals::operator*=(const double & zz) {
  if (comm_.rank() == 0) {
    aq_geovals_mult_f90(keyGeoVals_, zz);
  }
  return *this;
}
// -----------------------------------------------------------------------------
GeoVals & GeoVals::operator+=(const GeoVals & other) {
  if (comm_.rank() == 0) {
    aq_geovals_add_f90(keyGeoVals_, other.keyGeoVals_);
  }
  return *this;
}
// -----------------------------------------------------------------------------
GeoVals & GeoVals::operator-=(const GeoVals & other) {
  if (comm_.rank() == 0) {
    aq_geovals_diff_f90(keyGeoVals_, other.keyGeoVals_);
  }
  return *this;
}
// -----------------------------------------------------------------------------
GeoVals & GeoVals::operator*=(const GeoVals & other) {
  if (comm_.rank() == 0) {
    aq_geovals_schurmult_f90(keyGeoVals_, other.keyGeoVals_);
  }
  return *this;
}
// -----------------------------------------------------------------------------
double GeoVals::dot_product_with(const GeoVals & other) const {
  double zz;
  if (comm_.rank() == 0) {
    aq_geovals_dotprod_f90(keyGeoVals_, other.keyGeoVals_, zz);
  }
  comm_.broadcast(zz, 0);
  return zz;
}
// -----------------------------------------------------------------------------
void GeoVals::read(const eckit::Configuration & config) {
  if (comm_.rank() == 0) {
    aq_geovals_read_file_f90(keyGeoVals_, config);
  }
}
// -----------------------------------------------------------------------------
void GeoVals::write(const eckit::Configuration & config) const {
  if (comm_.rank() == 0) {
    aq_geovals_write_file_f90(keyGeoVals_, config);
  }
}
// -----------------------------------------------------------------------------
void GeoVals::print(std::ostream & os) const {
  if (comm_.rank() == 0) {
    int nobs;
    double zmin, zmax, zrms;
    aq_geovals_stats_f90(keyGeoVals_, nobs, zmin, zmax, zrms);
    std::ios_base::fmtflags f(os.flags());
    os << " nobs= " << nobs << std::scientific << std::setprecision(4)
       << "  Min=" << std::setw(12) << zmin
       << ", Max=" << std::setw(12) << zmax
       << ", RMS=" << std::setw(12) << zrms;
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
}
// -----------------------------------------------------------------------------
}  // namespace aq
