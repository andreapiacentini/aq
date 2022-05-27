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

#ifndef AQ_STATE_H_
#define AQ_STATE_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "aq/Fields.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace aq {
  class GeoVals;
  class Locations;
  class Geometry;
  class Increment;

/// AQ model state
// -----------------------------------------------------------------------------
class State : public util::Printable,
                private util::ObjectCounter<State> {
 public:
  static const std::string classname() {return "aq::State";}

/// Constructor, destructor
  State(const Geometry &, const oops::Variables &, const util::DateTime &);  // Is it used?
  State(const Geometry &, const eckit::Configuration &);
  State(const Geometry &, const State &);
  State(const State &);
  virtual ~State();
  State & operator=(const State &);

/// Interpolate full fields
  void changeResolution(const State & xx);

/// Interactions with Increment
  State & operator+=(const Increment &);

/// I/O and diagnostics
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const {return fields_->norm();}
  const util::DateTime & validTime() const {return fields_->time();}
  util::DateTime & validTime() {return fields_->time();}
  void updateTime(const util::Duration & dt) {fields_->updateTime(dt);}

/// Access to fields
  Fields & fields() {return *fields_;}
  const Fields & fields() const {return *fields_;}
  std::shared_ptr<const Geometry> geometry() const {
    return fields_->geometry();
  }
  const oops::Variables & variables() const {return fields_->variables();}

/// Serialization
  size_t serialSize() const;
  void serialize(std::vector<double> &) const;
  void deserialize(const std::vector<double> &, size_t &);

/// Other
  void zero();
  void accumul(const double &, const State &);

 private:
  void print(std::ostream &) const;
  std::unique_ptr<Fields> fields_;
};
// -----------------------------------------------------------------------------

}  // namespace aq

#endif  // AQ_STATE_H_
