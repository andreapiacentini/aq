/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2017-2019 UCAR.
 * (C) Copyright 2021-2022 CERFACS.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "aq/Increment.h"

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

#include "atlas/field.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

#include "aq/Covariance.h"
#include "aq/Fields.h"
#include "aq/Geometry.h"
#include "aq/GeoVals.h"
#include "aq/State.h"

namespace aq {

// -----------------------------------------------------------------------------
/// Constructor, destructor
// -----------------------------------------------------------------------------
Increment::Increment(const Geometry & resol, const oops::Variables & vars,
                         const util::DateTime & vt)
  : fields_(new Fields(resol, vars, vt)), vars_(vars)
{
  oops::Log::trace() << "Increment constructed." << std::endl;
}
// -----------------------------------------------------------------------------
Increment::Increment(const Geometry & resol, const Increment & other)
  : fields_(new Fields(*other.fields_, resol)), vars_(other.vars_)
{
  oops::Log::trace() << "Increment constructed from other." << std::endl;
}
// -----------------------------------------------------------------------------
Increment::Increment(const Increment & other, const bool copy)
  : fields_(new Fields(*other.fields_, copy)), vars_(other.vars_)
{
  oops::Log::trace() << "Increment copy-created." << std::endl;
}
// -----------------------------------------------------------------------------
Increment::Increment(const Increment & other)
  : fields_(new Fields(*other.fields_)), vars_(other.vars_)
{
  oops::Log::trace() << "Increment copy-created." << std::endl;
}
// -----------------------------------------------------------------------------
Increment::~Increment() {
  oops::Log::trace() << "Increment destructed" << std::endl;
}
// -----------------------------------------------------------------------------
/// Basic operators
// -----------------------------------------------------------------------------
void Increment::diff(const State & x1, const State & x2) {
  ASSERT(this->validTime() == x1.validTime());
  ASSERT(this->validTime() == x2.validTime());
  fields_->diff(x1.fields(), x2.fields());
}
// -----------------------------------------------------------------------------
Increment & Increment::operator=(const Increment & rhs) {
  *fields_ = *rhs.fields_;
  return *this;
}
// -----------------------------------------------------------------------------
Increment & Increment::operator+=(const Increment & dx) {
  ASSERT(this->validTime() == dx.validTime());
  *fields_ += *dx.fields_;
  return *this;
}
// -----------------------------------------------------------------------------
Increment & Increment::operator-=(const Increment & dx) {
  ASSERT(this->validTime() == dx.validTime());
  *fields_ -= *dx.fields_;
  return *this;
}
// -----------------------------------------------------------------------------
Increment & Increment::operator*=(const double & zz) {
  *fields_ *= zz;
  return *this;
}
// -----------------------------------------------------------------------------
void Increment::zero() {
  fields_->zero();
}
// -----------------------------------------------------------------------------
void Increment::zero(const util::DateTime & vt) {
  fields_->zero(vt);
}
// -----------------------------------------------------------------------------
void Increment::ones() {
  fields_->ones();
}
// -----------------------------------------------------------------------------
void Increment::axpy(const double & zz, const Increment & dx,
                       const bool check) {
  ASSERT(!check || this->validTime() == dx.validTime());
  fields_->axpy(zz, *dx.fields_);
}
// -----------------------------------------------------------------------------
void Increment::accumul(const double & zz, const State & xx) {
  fields_->axpy(zz, xx.fields());
}
// -----------------------------------------------------------------------------
void Increment::schur_product_with(const Increment & dx) {
  fields_->schur_product_with(*dx.fields_);
}
// -----------------------------------------------------------------------------
double Increment::dot_product_with(const Increment & other) const {
  return dot_product(*fields_, *other.fields_);
}
// -----------------------------------------------------------------------------
void Increment::random() {
  fields_->random();
}
// -----------------------------------------------------------------------------
void Increment::dirac(const eckit::Configuration & config) {
  fields_->zero();
  if (config.has("date")) {
    util::DateTime dd(config.getString("date"));
    if (this->validTime() == dd) fields_->dirac(config);
  } else if (config.has("dates")) {
    std::vector<std::string> dates(config.getStringVector("dates"));
    for (const auto date : dates) {
      util::DateTime dd(date);
      if (this->validTime() == dd) fields_->dirac(config);
    }
  } else {
    ABORT("missing date or dates in dirac configuration");
  }
}
// -----------------------------------------------------------------------------
/// ATLAS FieldSet
// -----------------------------------------------------------------------------
void Increment::toFieldSet(atlas::FieldSet & fset) const {
  fields_->toFieldSet(fset);
}
// -----------------------------------------------------------------------------
void Increment::toFieldSetAD(const atlas::FieldSet & fset) {
  fields_->toFieldSetAD(fset);
}
// -----------------------------------------------------------------------------
void Increment::fromFieldSet(const atlas::FieldSet & fset) {
  fields_->fromFieldSet(fset);
}
// -----------------------------------------------------------------------------
/// I/O and diagnostics
// -----------------------------------------------------------------------------
void Increment::read(const eckit::Configuration & files) {
  fields_->read(files);
}
// -----------------------------------------------------------------------------
void Increment::write(const eckit::Configuration & files) const {
  fields_->write(files);
}
// -----------------------------------------------------------------------------
std::vector<double> Increment::rmsByLevel(const std::string & var) const {
  return fields_->rmsByLevel(var);
}
// -----------------------------------------------------------------------------
/// Serialization
// -----------------------------------------------------------------------------
size_t Increment::serialSize() const {
  size_t nn = fields_->serialSize();
  return nn;
}
// -----------------------------------------------------------------------------
void Increment::serialize(std::vector<double> & vect) const {
  fields_->serialize(vect);
}
// -----------------------------------------------------------------------------
void Increment::deserialize(const std::vector<double> & vect, size_t & index) {
  fields_->deserialize(vect, index);
}
// -----------------------------------------------------------------------------
void Increment::print(std::ostream & os) const {
  os << std::endl << "  Valid time: " << validTime();
  os << *fields_;
}
// -----------------------------------------------------------------------------
oops::LocalIncrement Increment::getLocal(const GeometryIterator & iter) const {
  return fields_->getLocal(iter);
}
// -----------------------------------------------------------------------------
void Increment::setLocal(const oops::LocalIncrement & values, const GeometryIterator & iter) {
  fields_->setLocal(values, iter);
}
// -----------------------------------------------------------------------------

}  // namespace aq
