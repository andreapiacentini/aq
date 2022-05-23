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
