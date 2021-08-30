/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "aq/GomAQ.h"

#include <iomanip>

#include "aq/AqFortran.h"
#include "aq/LocationsAQ.h"
#include "aq/ObsSpaceAQ.h"
#include "eckit/config/Configuration.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

namespace aq {

// -----------------------------------------------------------------------------
GomAQ::GomAQ(const LocationsAQ & locs, const oops::Variables & vars,
             const std::vector<size_t> & sizes):
  vars_(vars)
{
  // gom_setup just creates and allocates the GeoVaLs object without filling
  // in values
  aq_gom_setup_f90(keyGom_, locs, vars_);
}
// -----------------------------------------------------------------------------
/*! AQ GeoVaLs Constructor with Config */

  GomAQ::GomAQ(const eckit::Configuration & config,
               const ObsSpaceAQ & ospace, const oops::Variables & vars):
  vars_(vars)
{
  aq_gom_create_f90(keyGom_, vars_);
  aq_gom_read_file_f90(keyGom_, config);
}
// -----------------------------------------------------------------------------
// Copy constructor
GomAQ::GomAQ(const GomAQ & other):
  vars_(other.vars_)
{
  aq_gom_create_f90(keyGom_, vars_);
  aq_gom_copy_f90(keyGom_, other.keyGom_);
}
// -----------------------------------------------------------------------------
GomAQ::~GomAQ() {
  aq_gom_delete_f90(keyGom_);
}
// -----------------------------------------------------------------------------
double GomAQ::rms() const {
  double zz;
  aq_gom_rms_f90(keyGom_, zz);
  return zz;
}
// -----------------------------------------------------------------------------
double GomAQ::normalizedrms(const GomAQ & rhs) const {
  GomAQ temp_GomAQ(*this);
  aq_gom_divide_f90(temp_GomAQ.keyGom_, rhs.keyGom_);
  return temp_GomAQ.rms();
}
// -----------------------------------------------------------------------------
void GomAQ::zero() {
  aq_gom_zero_f90(keyGom_);
}
// -----------------------------------------------------------------------------
void GomAQ::random() {
  aq_gom_random_f90(keyGom_);
}
// -----------------------------------------------------------------------------
GomAQ & GomAQ::operator=(const GomAQ & rhs) {
  const int keyGomRhs = rhs.keyGom_;
  aq_gom_copy_f90(keyGom_, keyGomRhs);
  return *this;
}
// -----------------------------------------------------------------------------
GomAQ & GomAQ::operator*=(const double & zz) {
  aq_gom_mult_f90(keyGom_, zz);
  return *this;
}
// -----------------------------------------------------------------------------
GomAQ & GomAQ::operator+=(const GomAQ & other) {
  aq_gom_add_f90(keyGom_, other.keyGom_);
  return *this;
}
// -----------------------------------------------------------------------------
GomAQ & GomAQ::operator-=(const GomAQ & other) {
  aq_gom_diff_f90(keyGom_, other.keyGom_);
  return *this;
}
// -----------------------------------------------------------------------------
GomAQ & GomAQ::operator*=(const GomAQ & other) {
  aq_gom_schurmult_f90(keyGom_, other.keyGom_);
  return *this;
}
// -----------------------------------------------------------------------------
double GomAQ::dot_product_with(const GomAQ & other) const {
  double zz;
  aq_gom_dotprod_f90(keyGom_, other.keyGom_, zz);
  return zz;
}
// -----------------------------------------------------------------------------
void GomAQ::read(const eckit::Configuration & config) {
  aq_gom_read_file_f90(keyGom_, config);
}
// -----------------------------------------------------------------------------
void GomAQ::write(const eckit::Configuration & config) const {
  aq_gom_write_file_f90(keyGom_, config);
}
// -----------------------------------------------------------------------------
void GomAQ::print(std::ostream & os) const {
  int nobs;
  double zmin, zmax, zrms;
  aq_gom_stats_f90(keyGom_, nobs, zmin, zmax, zrms);
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

    aq_gom_maxloc_f90(keyGom_, mxval, iloc, maxvar);

    oops::Log::debug() << "GomAQ: Maximum Value = " << std::setprecision(4)
                       << mxval << " at location = " << iloc
                       << " and variable = " << maxvar << std::endl;
  }
}
// -----------------------------------------------------------------------------
}  // namespace aq
