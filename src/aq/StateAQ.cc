/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "aq/StateAQ.h"

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

#include "aq/FieldsAQ.h"
#include "aq/GeometryAQ.h"
#include "aq/GomAQ.h"
#include "aq/IncrementAQ.h"
#include "aq/LocationsAQ.h"


namespace aq {

// -----------------------------------------------------------------------------
/// Constructor, destructor
// -----------------------------------------------------------------------------
StateAQ::StateAQ(const GeometryAQ & resol, const oops::Variables & vars,
                 const util::DateTime & vt)
  : fields_(new FieldsAQ(resol, vars, vt))
{
  oops::Log::trace() << "StateAQ::StateAQ created." << std::endl;
}
// -----------------------------------------------------------------------------
StateAQ::StateAQ(const GeometryAQ & resol, const eckit::Configuration & file)
  : fields_()
{
  oops::Variables vars(file, "state variables");
  oops::Log::trace() << "StateAQ::StateAQ variables: " << vars << std::endl;
  fields_.reset(new FieldsAQ(resol, vars, util::DateTime(file.getString("date"))));
  if (file.has("analytic_init")) {
    fields_->analytic_init(file);
  } else if (file.has("read_from_file")) {
    const int read_from_file = file.getInt("read_from_file");
    if (read_from_file == 0) {
       fields_->analytic_init(file);
    } else if (read_from_file == 1) {
      fields_->read(file);
    }
  } else {
    fields_->read(file);
  }

  ASSERT(fields_);
  oops::Log::trace() << "StateAQ::StateAQ created and read in." << std::endl;
}
// -----------------------------------------------------------------------------
StateAQ::StateAQ(const GeometryAQ & resol, const StateAQ & other)
  : fields_(new FieldsAQ(*other.fields_, resol))
{
  ASSERT(fields_);
  oops::Log::trace() << "StateAQ::StateAQ created by interpolation." << std::endl;
}
// -----------------------------------------------------------------------------
StateAQ::StateAQ(const StateAQ & other)
  : fields_(new FieldsAQ(*other.fields_))
{
  ASSERT(fields_);
  oops::Log::trace() << "StateAQ::StateAQ copied." << std::endl;
}
// -----------------------------------------------------------------------------
StateAQ::~StateAQ() {
  oops::Log::trace() << "StateAQ::StateAQ destructed." << std::endl;
}
// -----------------------------------------------------------------------------
/// Basic operators
// -----------------------------------------------------------------------------
StateAQ & StateAQ::operator=(const StateAQ & rhs) {
  ASSERT(fields_);
  *fields_ = *rhs.fields_;
  return *this;
}
// -----------------------------------------------------------------------------
/// Interpolate full fields
// -----------------------------------------------------------------------------
void StateAQ::changeResolution(const StateAQ & other) {
  fields_->changeResolution(*other.fields_);
  oops::Log::trace() << "StateAQ interpolated" << std::endl;
}
// -----------------------------------------------------------------------------
/// Interactions with Increments
// -----------------------------------------------------------------------------
StateAQ & StateAQ::operator+=(const IncrementAQ & dx) {
  ASSERT(this->validTime() == dx.validTime());
  ASSERT(fields_);
  fields_->add(dx.fields());
  return *this;
}
// -----------------------------------------------------------------------------
/// I/O and diagnostics
// -----------------------------------------------------------------------------
void StateAQ::read(const eckit::Configuration & files) {
  fields_->read(files);
}
// -----------------------------------------------------------------------------
void StateAQ::write(const eckit::Configuration & files) const {
  fields_->write(files);
}
// -----------------------------------------------------------------------------
/// Serialization
// -----------------------------------------------------------------------------
size_t StateAQ::serialSize() const {
  size_t nn = fields_->serialSize();
  return nn;
}
// -----------------------------------------------------------------------------
void StateAQ::serialize(std::vector<double> & vect) const {
  fields_->serialize(vect);
}
// -----------------------------------------------------------------------------
void StateAQ::deserialize(const std::vector<double> & vect, size_t & index) {
  fields_->deserialize(vect, index);
}
// -----------------------------------------------------------------------------
void StateAQ::print(std::ostream & os) const {
  os << std::endl << "  Valid time: " << validTime();
  os << *fields_;
}
// -----------------------------------------------------------------------------
/// For accumulator
// -----------------------------------------------------------------------------
void StateAQ::zero() {
  fields_->zero();
}
// -----------------------------------------------------------------------------
void StateAQ::accumul(const double & zz, const StateAQ & xx) {
  fields_->axpy(zz, *xx.fields_);
}

}  // namespace aq
