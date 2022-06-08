/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2021-2022 CERFACS.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "aq/Fields.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <map>
#include <numeric>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/utils/StringTools.h"

#include "aq/aq_fields_interface.h"
#include "aq/Geometry.h"
#include "aq/GeoVals.h"
#include "aq/Locations.h"

#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

#include "oops/util/abor1_cpp.h"

// -----------------------------------------------------------------------------
namespace aq {
// -----------------------------------------------------------------------------
Fields::Fields(const Geometry & geom, const oops::Variables & vars,
                  const util::DateTime & time):
  geom_(geom), vars_(vars), time_(time)
{
  aq_fields_create_f90(keyFlds_, geom_.toFortran(), vars_);
}
// -----------------------------------------------------------------------------
Fields::Fields(const Fields & other, const bool copy)
  : geom_(other.geom_), vars_(other.vars_), time_(other.time_)
{
  aq_fields_create_from_other_f90(keyFlds_, other.keyFlds_);
  if (copy) {
    aq_fields_copy_f90(keyFlds_, other.keyFlds_);
  }
}
// -----------------------------------------------------------------------------
Fields::Fields(const Fields & other)
  : geom_(other.geom_), vars_(other.vars_), time_(other.time_)
{
  aq_fields_create_from_other_f90(keyFlds_, other.keyFlds_);
  aq_fields_copy_f90(keyFlds_, other.keyFlds_);
}
// -----------------------------------------------------------------------------
Fields::Fields(const Fields & other, const Geometry & geom)
  : geom_(geom), vars_(other.vars_), time_(other.time_)
{
  int world_rank = oops::mpi::world().rank();
  if (world_rank == 0) {
    oops::Log::debug() << " Change of resolution not yet implemented. Copy instead. " << std::endl;
  }
  // TODO(Emanuele and Andrea): Implement the change of resolution.
  // aq_fields_create_f90(keyFlds_, geom_.toFortran(), vars_);
  // aq_fields_change_resol_f90(keyFlds_, other.keyFlds_);
  aq_fields_create_from_other_f90(keyFlds_, other.keyFlds_);
  aq_fields_copy_f90(keyFlds_, other.keyFlds_);
}
// -----------------------------------------------------------------------------
Fields::Fields(const Fields & other, const oops::Variables & vars)
  : geom_(other.geom_), vars_(vars), time_(other.time_)
{
// TODO(Benjamin): delete that ?
  aq_fields_create_f90(keyFlds_, geom_.toFortran(), vars_);
  aq_fields_copy_f90(keyFlds_, other.keyFlds_);
}
// -----------------------------------------------------------------------------
Fields::~Fields() {
  aq_fields_delete_f90(keyFlds_);
}
// -----------------------------------------------------------------------------
Fields & Fields::operator=(const Fields & rhs) {
  aq_fields_copy_f90(keyFlds_, rhs.keyFlds_);
  time_ = rhs.time_;
  return *this;
}
// -----------------------------------------------------------------------------
Fields & Fields::operator+=(const Fields & rhs) {
  aq_fields_self_add_f90(keyFlds_, rhs.keyFlds_);
  return *this;
}
// -----------------------------------------------------------------------------
Fields & Fields::operator-=(const Fields & rhs) {
  aq_fields_self_sub_f90(keyFlds_, rhs.keyFlds_);
  return *this;
}
// -----------------------------------------------------------------------------
Fields & Fields::operator*=(const double & zz) {
  aq_fields_self_mul_f90(keyFlds_, zz);
  return *this;
}
// -----------------------------------------------------------------------------
void Fields::zero() {
  aq_fields_zero_f90(keyFlds_);
}
// -----------------------------------------------------------------------------
void Fields::zero(const util::DateTime & time) {
  aq_fields_zero_f90(keyFlds_);
  time_ = time;
}
// -----------------------------------------------------------------------------
void Fields::ones() {
  aq_fields_ones_f90(keyFlds_);
}
// -----------------------------------------------------------------------------
void Fields::axpy(const double & zz, const Fields & rhs) {
  aq_fields_axpy_f90(keyFlds_, zz, rhs.keyFlds_);
}
// -----------------------------------------------------------------------------
double Fields::dot_product_with(const Fields & fld2) const {
  double zz;
  aq_fields_dot_prod_f90(keyFlds_, fld2.keyFlds_, zz);
  return zz;
}
// -----------------------------------------------------------------------------
void Fields::schur_product_with(const Fields & dx) {
    aq_fields_self_schur_f90(keyFlds_, dx.keyFlds_);
}
// -----------------------------------------------------------------------------
void Fields::random() {
  aq_fields_random_f90(keyFlds_);
}
// -----------------------------------------------------------------------------
void Fields::dirac(const eckit::Configuration & config) {
  aq_fields_dirac_f90(keyFlds_, config);
}
// -----------------------------------------------------------------------------
void Fields::changeResolution(const Fields & other) {
  // aq_fields_change_resol_f90(keyFlds_, other.keyFlds_);
}
// -----------------------------------------------------------------------------
void Fields::add(const Fields & rhs) {
  // TODO(Emanuele and Andrea): Implement change of resolution
  // Fields rhs_myres(rhs, geom_);
  // aq_fields_add_incr_f90(keyFlds_, rhs_myres.keyFlds_);
  aq_fields_add_incr_f90(keyFlds_, rhs.keyFlds_);
}
// -----------------------------------------------------------------------------
void Fields::diff(const Fields & x1, const Fields & x2) {
  // TODO(Emanuele and Andrea): Implement change of resolution
  // Fields x1_myres(x1, geom_);
  // Fields x2_myres(x2, geom_);
  // aq_fields_diff_incr_f90(keyFlds_, x1_myres.keyFlds_, x2_myres.keyFlds_);
  aq_fields_diff_incr_f90(keyFlds_, x1.keyFlds_, x2.keyFlds_);
}
// -----------------------------------------------------------------------------
void Fields::setAtlas(atlas::FieldSet * afieldset) const {
  aq_fields_set_atlas_f90(keyFlds_, vars_, afieldset->get());
}
// -----------------------------------------------------------------------------
void Fields::toAtlas(atlas::FieldSet * afieldset) const {
  aq_fields_to_atlas_f90(keyFlds_, vars_, afieldset->get());
}
// -----------------------------------------------------------------------------
void Fields::fromAtlas(atlas::FieldSet * afieldset) {
  oops::Log::debug() << "fromAtlas is a passive empty call for Fields coded as atlas fieldsets"
                     << std::endl;
  // aq_fields_from_atlas_f90(keyFlds_, vars_, afieldset->get());
}
// -----------------------------------------------------------------------------
void Fields::read(const eckit::Configuration & config) {
  aq_fields_read_file_f90(keyFlds_, config, time_);
}
// -----------------------------------------------------------------------------
void Fields::analytic_init(const eckit::Configuration & config) {
  aq_fields_analytic_init_f90(keyFlds_, config);
}
// -----------------------------------------------------------------------------
void Fields::write(const eckit::Configuration & config) const {
  aq_fields_write_file_f90(keyFlds_, config, time_);
}
// -----------------------------------------------------------------------------
double Fields::norm() const {
  double zz = 0.0;
  aq_fields_rms_f90(keyFlds_, zz);
  return zz;
}
// -----------------------------------------------------------------------------
void Fields::print(std::ostream & os) const {
  eckit::LocalConfiguration flds_info_;
  aq_fields_info_f90(keyFlds_, flds_info_);
  std::string name = flds_info_.getString("name");
  int nx = flds_info_.getInt("geometry.nx");
  int ny = flds_info_.getInt("geometry.ny");
  int nz = flds_info_.getInt("geometry.nz");
  int mod_levels = flds_info_.getInt("geometry.modlev");
  std::string orientation = flds_info_.getString("geometry.orientation");
  std::string domname = flds_info_.getString("geometry.domname");
  std::string model = flds_info_.getString("geometry.model");
  os << std::endl << "  Fields " << name << " for model " << model << std::endl;
  os << "  on domain " << domname << " (nx=" << nx << ", ny=" << ny <<
    ", lev=" << mod_levels << ")" << std::endl;
  os << "  " << nz << " vertical levels (orientation " << orientation << ")" << std::endl;
  std::vector<std::string> species = flds_info_.getStringVector("state variables");
  os << "  with " << species.size() << " variables" << std::endl;
  double min, max, mean, stddev;
  for (int i = 0; i < static_cast<int>(species.size()); i++) {
    min = flds_info_.getDouble("statistics."+eckit::StringTools::trim(species[i])+".min");
    max = flds_info_.getDouble("statistics."+eckit::StringTools::trim(species[i])+".max");
    mean = flds_info_.getDouble("statistics."+eckit::StringTools::trim(species[i])+".mean");
    stddev = flds_info_.getDouble("statistics."+eckit::StringTools::trim(species[i])+".stddev");
    os << "  Var. " << eckit::StringTools::trim(species[i])
       << ": min = " << min << "; max = " << max
       << "; mean = " << mean << "; stddev = " << stddev << std::endl;
  }
}
// -----------------------------------------------------------------------------
size_t Fields::serialSize() const {
  int nnInt;
  aq_fields_serialsize_f90(keyFlds_, nnInt);
  size_t nn = nnInt;
  nn += time_.serialSize();
  return nn;
}
// -----------------------------------------------------------------------------
void Fields::serialize(std::vector<double> & vect)  const {
  size_t size_fld = this->serialSize() - time_.serialSize();

  // Allocate space for fld, xb and qb
  std::vector<double> v_fld(size_fld, 0);

  // Serialize the field
  aq_fields_serialize_f90(keyFlds_, static_cast<int>(size_fld), v_fld.data());
  vect.insert(vect.end(), v_fld.begin(), v_fld.end());

  // Serialize the date and time
  time_.serialize(vect);
}
// -----------------------------------------------------------------------------
void Fields::deserialize(const std::vector<double> & vect, size_t & index) {
  int indexInt = static_cast<int>(index);
  aq_fields_deserialize_f90(keyFlds_, static_cast<int>(vect.size()), vect.data(), indexInt);
  index = static_cast<size_t>(indexInt);
  time_.deserialize(vect, index);
}
// -----------------------------------------------------------------------------
}  // namespace aq
