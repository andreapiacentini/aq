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

#include "aq/State.h"

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

#include "aq/Fields.h"
#include "aq/Geometry.h"
#include "aq/GeoVals.h"
#include "aq/Increment.h"
#include "aq/Locations.h"


namespace aq {

// -----------------------------------------------------------------------------
/// Constructor, destructor
// -----------------------------------------------------------------------------
State::State(const Geometry & resol, const oops::Variables & vars,
                 const util::DateTime & vt)
  : fields_(new Fields(resol, vars, vt))
{
  oops::Log::trace() << "State::State created." << std::endl;
}
// -----------------------------------------------------------------------------
State::State(const Geometry & resol, const eckit::Configuration & file)
  : fields_()
{
  oops::Variables vars(file, "state variables");
  oops::Log::trace() << "State::State variables: " << vars << std::endl;
  fields_.reset(new Fields(resol, vars, util::DateTime(file.getString("date"))));
  if (file.has("analytic init")) {
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
  oops::Log::trace() << "State::State created and read in." << std::endl;
}
// -----------------------------------------------------------------------------
State::State(const Geometry & resol, const State & other)
  : fields_(new Fields(*other.fields_, resol))
{
  ASSERT(fields_);
  oops::Log::trace() << "State::State created by interpolation." << std::endl;
}
// -----------------------------------------------------------------------------
State::State(const State & other)
  : fields_(new Fields(*other.fields_))
{
  ASSERT(fields_);
  oops::Log::trace() << "State::State copied." << std::endl;
}
// -----------------------------------------------------------------------------
State::~State() {
  oops::Log::trace() << "State::State destructed." << std::endl;
}
// -----------------------------------------------------------------------------
/// Basic operators
// -----------------------------------------------------------------------------
State & State::operator=(const State & rhs) {
  ASSERT(fields_);
  *fields_ = *rhs.fields_;
  return *this;
}
// -----------------------------------------------------------------------------
/// Interpolate full fields
// -----------------------------------------------------------------------------
void State::changeResolution(const State & other) {
  fields_->changeResolution(*other.fields_);
  oops::Log::trace() << "State interpolated" << std::endl;
}
// -----------------------------------------------------------------------------
/// Interactions with Increments
// -----------------------------------------------------------------------------
State & State::operator+=(const Increment & dx) {
  ASSERT(this->validTime() == dx.validTime());
  ASSERT(fields_);
  fields_->add(dx.fields());
  return *this;
}
// -----------------------------------------------------------------------------
/// I/O and diagnostics
// -----------------------------------------------------------------------------
void State::read(const eckit::Configuration & files) {
  fields_->read(files);
}
// -----------------------------------------------------------------------------
void State::write(const eckit::Configuration & files) const {
  fields_->write(files);
}
// -----------------------------------------------------------------------------
/// Serialization
// -----------------------------------------------------------------------------
size_t State::serialSize() const {
  size_t nn = fields_->serialSize();
  return nn;
}
// -----------------------------------------------------------------------------
void State::serialize(std::vector<double> & vect) const {
  fields_->serialize(vect);
}
// -----------------------------------------------------------------------------
void State::deserialize(const std::vector<double> & vect, size_t & index) {
  fields_->deserialize(vect, index);
}
// -----------------------------------------------------------------------------
void State::print(std::ostream & os) const {
  os << std::endl << "  Valid time: " << validTime();
  os << *fields_;
}
// -----------------------------------------------------------------------------
/// ATLAS FieldSet
// -----------------------------------------------------------------------------
void State::toFieldSet(atlas::FieldSet & fset) const {
  fields_->toFieldSet(fset);
}
// -----------------------------------------------------------------------------
void State::fromFieldSet(const atlas::FieldSet & fset) {
  fields_->fromFieldSet(fset);
}
// -----------------------------------------------------------------------------
/// For accumulator
// -----------------------------------------------------------------------------
void State::zero() {
  fields_->zero();
}
// -----------------------------------------------------------------------------
void State::accumul(const double & zz, const State & xx) {
  fields_->axpy(zz, *xx.fields_);
}

}  // namespace aq
