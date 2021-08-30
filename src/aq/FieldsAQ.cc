/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "aq/FieldsAQ.h"

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

#include "aq/AqFortran.h"
#include "aq/GeometryAQ.h"
#include "aq/GomAQ.h"
#include "aq/LocationsAQ.h"

#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

#include "oops/util/abor1_cpp.h"

// -----------------------------------------------------------------------------
namespace aq {
// -----------------------------------------------------------------------------
FieldsAQ::FieldsAQ(const GeometryAQ & geom, const oops::Variables & vars,
                  const util::DateTime & time):
  geom_(new GeometryAQ(geom)), vars_(vars), time_(time)
{
  aq_fields_create_f90(keyFlds_, geom_->toFortran(), vars_);
}
// -----------------------------------------------------------------------------
FieldsAQ::FieldsAQ(const FieldsAQ & other, const bool copy)
  : vars_(other.vars_), time_(other.time_)
{
  aq_fields_create_from_other_f90(keyFlds_, other.keyFlds_);
  if (copy) {
    aq_fields_copy_f90(keyFlds_, other.keyFlds_);
  }
}
// -----------------------------------------------------------------------------
FieldsAQ::FieldsAQ(const FieldsAQ & other)
  : vars_(other.vars_), time_(other.time_)
{
  aq_fields_create_from_other_f90(keyFlds_, other.keyFlds_);
  aq_fields_copy_f90(keyFlds_, other.keyFlds_);
}
// -----------------------------------------------------------------------------
FieldsAQ::FieldsAQ(const FieldsAQ & other, const GeometryAQ & geom)
  : geom_(new GeometryAQ(geom)), vars_(other.vars_), time_(other.time_)
{
  int world_rank = oops::mpi::world().rank();
  if (world_rank == 0) {
    oops::Log::debug() << " Change of resolution not yet implemented. Copy instead. " << std::endl;
  }
  // TODO(Emanuele and Andrea): Implement the change of resolution.
  // aq_fields_create_f90(keyFlds_, geom_->toFortran(), vars_);
  // aq_fields_change_resol_f90(keyFlds_, other.keyFlds_);
  aq_fields_create_from_other_f90(keyFlds_, other.keyFlds_);
  aq_fields_copy_f90(keyFlds_, other.keyFlds_);
}
// -----------------------------------------------------------------------------
FieldsAQ::FieldsAQ(const FieldsAQ & other, const oops::Variables & vars)
  : geom_(other.geom_), vars_(vars), time_(other.time_)
{
// TODO(Benjamin): delete that ?
  aq_fields_create_f90(keyFlds_, geom_->toFortran(), vars_);
  aq_fields_copy_f90(keyFlds_, other.keyFlds_);
}
// -----------------------------------------------------------------------------
FieldsAQ::~FieldsAQ() {
  aq_fields_delete_f90(keyFlds_);
}
// -----------------------------------------------------------------------------
FieldsAQ & FieldsAQ::operator=(const FieldsAQ & rhs) {
  aq_fields_copy_f90(keyFlds_, rhs.keyFlds_);
  time_ = rhs.time_;
  return *this;
}
// -----------------------------------------------------------------------------
FieldsAQ & FieldsAQ::operator+=(const FieldsAQ & rhs) {
  aq_fields_self_add_f90(keyFlds_, rhs.keyFlds_);
  return *this;
}
// -----------------------------------------------------------------------------
FieldsAQ & FieldsAQ::operator-=(const FieldsAQ & rhs) {
  aq_fields_self_sub_f90(keyFlds_, rhs.keyFlds_);
  return *this;
}
// -----------------------------------------------------------------------------
FieldsAQ & FieldsAQ::operator*=(const double & zz) {
  aq_fields_self_mul_f90(keyFlds_, zz);
  return *this;
}
// -----------------------------------------------------------------------------
void FieldsAQ::zero() {
  aq_fields_zero_f90(keyFlds_);
}
// -----------------------------------------------------------------------------
void FieldsAQ::zero(const util::DateTime & time) {
  aq_fields_zero_f90(keyFlds_);
  time_ = time;
}
// -----------------------------------------------------------------------------
void FieldsAQ::ones() {
  aq_fields_ones_f90(keyFlds_);
}
// -----------------------------------------------------------------------------
void FieldsAQ::axpy(const double & zz, const FieldsAQ & rhs) {
  aq_fields_axpy_f90(keyFlds_, zz, rhs.keyFlds_);
}
// -----------------------------------------------------------------------------
double FieldsAQ::dot_product_with(const FieldsAQ & fld2) const {
  double zz;
  aq_fields_dot_prod_f90(keyFlds_, fld2.keyFlds_, zz);
  return zz;
}
// -----------------------------------------------------------------------------
void FieldsAQ::schur_product_with(const FieldsAQ & dx) {
    aq_fields_self_schur_f90(keyFlds_, dx.keyFlds_);
}
// -----------------------------------------------------------------------------
void FieldsAQ::random() {
  aq_fields_random_f90(keyFlds_, vars_);
}
// -----------------------------------------------------------------------------
void FieldsAQ::dirac(const eckit::Configuration & config) {
  aq_fields_dirac_f90(keyFlds_, config);
}
// -----------------------------------------------------------------------------
void FieldsAQ::changeResolution(const FieldsAQ & other) {
  // aq_fields_change_resol_f90(keyFlds_, other.keyFlds_);
}
// -----------------------------------------------------------------------------
void FieldsAQ::add(const FieldsAQ & rhs) {
  // TODO(Emanuele and Andrea): Implement change of resolution
  // FieldsAQ rhs_myres(rhs, *geom_);
  // aq_fields_add_incr_f90(keyFlds_, rhs_myres.keyFlds_);
  aq_fields_add_incr_f90(keyFlds_, rhs.keyFlds_);
}
// -----------------------------------------------------------------------------
void FieldsAQ::diff(const FieldsAQ & x1, const FieldsAQ & x2) {
  // TODO(Emanuele and Andrea): Implement change of resolution
  // FieldsAQ x1_myres(x1, *geom_);
  // FieldsAQ x2_myres(x2, *geom_);
  // aq_fields_diff_incr_f90(keyFlds_, x1_myres.keyFlds_, x2_myres.keyFlds_);
  aq_fields_diff_incr_f90(keyFlds_, x1.keyFlds_, x2.keyFlds_);
}
// -----------------------------------------------------------------------------
void FieldsAQ::setAtlas(atlas::FieldSet * afieldset) const {
  aq_fields_set_atlas_f90(keyFlds_, vars_, afieldset->get());
}
// -----------------------------------------------------------------------------
void FieldsAQ::toAtlas(atlas::FieldSet * afieldset) const {
  aq_fields_to_atlas_f90(keyFlds_, vars_, afieldset->get());
}
// -----------------------------------------------------------------------------
void FieldsAQ::fromAtlas(atlas::FieldSet * afieldset) {
  oops::Log::debug() << "fromAtlas is a passive empty call for Fields coded as atlas fieldsets"
                     << std::endl;
  // aq_fields_from_atlas_f90(keyFlds_, vars_, afieldset->get());
}
// -----------------------------------------------------------------------------
void FieldsAQ::read(const eckit::Configuration & config) {
  aq_fields_read_file_f90(keyFlds_, config, time_);
}
// -----------------------------------------------------------------------------
void FieldsAQ::analytic_init(const eckit::Configuration & config) {
  aq_fields_analytic_init_f90(keyFlds_, config);
}
// -----------------------------------------------------------------------------
void FieldsAQ::write(const eckit::Configuration & config) const {
  aq_fields_write_file_f90(keyFlds_, config, time_);
}
// -----------------------------------------------------------------------------
double FieldsAQ::norm() const {
  double zz = 0.0;
  aq_fields_rms_f90(keyFlds_, zz);
  return zz;
}
// -----------------------------------------------------------------------------
void FieldsAQ::print(std::ostream & os) const {
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
  for (int i = 0; i < species.size(); i++) {
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
oops::LocalIncrement FieldsAQ::getLocal(const GeometryAQIterator & iter) const {
  // int nx, ny, nz;
  // aq_fields_sizes_f90(keyFlds_, nx, ny, nz);
  // std::vector<int> varlens(vars_.size());
  // for (unsigned int ii = 0; ii < vars_.size(); ii++) {
  //   varlens[ii] = nz;
  // }
  // int lenvalues = std::accumulate(varlens.begin(), varlens.end(), 0);
  // std::vector<double> values(lenvalues);
  // aq_fields_getpoint_f90(keyFlds_, iter.toFortran(), values.size(), values[0]);
  // return oops::LocalIncrement(vars_, values, varlens);
}
// -----------------------------------------------------------------------------
void FieldsAQ::setLocal(const oops::LocalIncrement & x, const GeometryAQIterator & iter) {
  // const std::vector<double> vals = x.getVals();
  // aq_fields_setpoint_f90(keyFlds_, iter.toFortran(), vals.size(), vals[0]);
}
// -----------------------------------------------------------------------------
size_t FieldsAQ::serialSize() const {
  int nnInt;
  aq_fields_serialsize_f90(keyFlds_, nnInt);
  size_t nn = nnInt;
  nn += time_.serialSize();
  return nn;
}
// -----------------------------------------------------------------------------
void FieldsAQ::serialize(std::vector<double> & vect)  const {
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
void FieldsAQ::deserialize(const std::vector<double> & vect, size_t & index) {
  int indexInt = static_cast<int>(index);
  aq_fields_deserialize_f90(keyFlds_, static_cast<int>(vect.size()), vect.data(), indexInt);
  index = static_cast<size_t>(indexInt);
  time_.deserialize(vect, index);
}
// -----------------------------------------------------------------------------
}  // namespace aq
