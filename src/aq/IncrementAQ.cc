/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2017-2019 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "model/IncrementAQ.h"

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

#include "model/ErrorCovarianceAQ.h"
#include "model/FieldsAQ.h"
#include "model/GeometryAQ.h"
#include "model/GomAQ.h"
#include "model/LocationsAQ.h"
#include "model/StateAQ.h"

namespace aq {

// -----------------------------------------------------------------------------
/// Constructor, destructor
// -----------------------------------------------------------------------------
IncrementAQ::IncrementAQ(const GeometryAQ & resol, const oops::Variables & vars,
                         const util::DateTime & vt)
  : fields_(new FieldsAQ(resol, vars, vt))
{
  fields_->zero();
  oops::Log::trace() << "IncrementAQ constructed." << std::endl;
}
// -----------------------------------------------------------------------------
IncrementAQ::IncrementAQ(const GeometryAQ & resol, const IncrementAQ & other)
  : fields_(new FieldsAQ(*other.fields_, resol))
{
  oops::Log::trace() << "IncrementAQ constructed from other." << std::endl;
}
// -----------------------------------------------------------------------------
IncrementAQ::IncrementAQ(const IncrementAQ & other, const bool copy)
  : fields_(new FieldsAQ(*other.fields_, copy))
{
  oops::Log::trace() << "IncrementAQ copy-created." << std::endl;
}
// -----------------------------------------------------------------------------
IncrementAQ::IncrementAQ(const IncrementAQ & other)
  : fields_(new FieldsAQ(*other.fields_))
{
  oops::Log::trace() << "IncrementAQ copy-created." << std::endl;
}
// -----------------------------------------------------------------------------
IncrementAQ::~IncrementAQ() {
  oops::Log::trace() << "IncrementAQ destructed" << std::endl;
}
// -----------------------------------------------------------------------------
/// Basic operators
// -----------------------------------------------------------------------------
void IncrementAQ::diff(const StateAQ & x1, const StateAQ & x2) {
  ASSERT(this->validTime() == x1.validTime());
  ASSERT(this->validTime() == x2.validTime());
  fields_->diff(x1.fields(), x2.fields());
}
// -----------------------------------------------------------------------------
IncrementAQ & IncrementAQ::operator=(const IncrementAQ & rhs) {
  *fields_ = *rhs.fields_;
  return *this;
}
// -----------------------------------------------------------------------------
IncrementAQ & IncrementAQ::operator+=(const IncrementAQ & dx) {
  ASSERT(this->validTime() == dx.validTime());
  *fields_ += *dx.fields_;
  return *this;
}
// -----------------------------------------------------------------------------
IncrementAQ & IncrementAQ::operator-=(const IncrementAQ & dx) {
  ASSERT(this->validTime() == dx.validTime());
  *fields_ -= *dx.fields_;
  return *this;
}
// -----------------------------------------------------------------------------
IncrementAQ & IncrementAQ::operator*=(const double & zz) {
  *fields_ *= zz;
  return *this;
}
// -----------------------------------------------------------------------------
void IncrementAQ::zero() {
  fields_->zero();
}
// -----------------------------------------------------------------------------
void IncrementAQ::zero(const util::DateTime & vt) {
  fields_->zero(vt);
}
// -----------------------------------------------------------------------------
void IncrementAQ::ones() {
  fields_->ones();
}
// -----------------------------------------------------------------------------
void IncrementAQ::axpy(const double & zz, const IncrementAQ & dx,
                       const bool check) {
  ASSERT(!check || this->validTime() == dx.validTime());
  fields_->axpy(zz, *dx.fields_);
}
// -----------------------------------------------------------------------------
void IncrementAQ::accumul(const double & zz, const StateAQ & xx) {
  fields_->axpy(zz, xx.fields());
}
// -----------------------------------------------------------------------------
void IncrementAQ::schur_product_with(const IncrementAQ & dx) {
  fields_->schur_product_with(*dx.fields_);
}
// -----------------------------------------------------------------------------
double IncrementAQ::dot_product_with(const IncrementAQ & other) const {
  return dot_product(*fields_, *other.fields_);
}
// -----------------------------------------------------------------------------
void IncrementAQ::random() {
  fields_->random();
}
// -----------------------------------------------------------------------------
void IncrementAQ::dirac(const eckit::Configuration & config) {
  fields_->zero();
  util::DateTime dd(config.getString("date"));
  if (this->validTime() == dd) fields_->dirac(config);
}
// -----------------------------------------------------------------------------
/// ATLAS FieldSet
// -----------------------------------------------------------------------------
void IncrementAQ::setAtlas(atlas::FieldSet * afieldset) const {
  fields_->setAtlas(afieldset);
}
// -----------------------------------------------------------------------------
void IncrementAQ::toAtlas(atlas::FieldSet * afieldset) const {
  fields_->toAtlas(afieldset);
}
// -----------------------------------------------------------------------------
void IncrementAQ::fromAtlas(atlas::FieldSet * afieldset) {
  fields_->fromAtlas(afieldset);
}
// -----------------------------------------------------------------------------
/// I/O and diagnostics
// -----------------------------------------------------------------------------
void IncrementAQ::read(const eckit::Configuration & files) {
  fields_->read(files);
}
// -----------------------------------------------------------------------------
void IncrementAQ::write(const eckit::Configuration & files) const {
  fields_->write(files);
}
// -----------------------------------------------------------------------------
/// Serialization
// -----------------------------------------------------------------------------
size_t IncrementAQ::serialSize() const {
  size_t nn = fields_->serialSize();
  return nn;
}
// -----------------------------------------------------------------------------
void IncrementAQ::serialize(std::vector<double> & vect) const {
  fields_->serialize(vect);
}
// -----------------------------------------------------------------------------
void IncrementAQ::deserialize(const std::vector<double> & vect, size_t & index) {
  fields_->deserialize(vect, index);
}
// -----------------------------------------------------------------------------
void IncrementAQ::print(std::ostream & os) const {
  os << std::endl << "  Valid time: " << validTime();
  os << *fields_;
}
// -----------------------------------------------------------------------------
oops::LocalIncrement IncrementAQ::getLocal(const GeometryAQIterator & iter) const {
  return fields_->getLocal(iter);
}
// -----------------------------------------------------------------------------
void IncrementAQ::setLocal(const oops::LocalIncrement & values, const GeometryAQIterator & iter) {
  fields_->setLocal(values, iter);
}
// -----------------------------------------------------------------------------

}  // namespace aq
