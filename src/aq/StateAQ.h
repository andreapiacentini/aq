/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef AQ_STATEAQ_H_
#define AQ_STATEAQ_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "aq/FieldsAQ.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace aq {
  class GomAQ;
  class LocationsAQ;
  class GeometryAQ;
  class IncrementAQ;

/// AQ model state
// -----------------------------------------------------------------------------
class StateAQ : public util::Printable,
                private util::ObjectCounter<StateAQ> {
 public:
  static const std::string classname() {return "aq::StateAQ";}

/// Constructor, destructor
  StateAQ(const GeometryAQ &, const oops::Variables &, const util::DateTime &);  // Is it used?
  StateAQ(const GeometryAQ &, const eckit::Configuration &);
  StateAQ(const GeometryAQ &, const StateAQ &);
  StateAQ(const StateAQ &);
  virtual ~StateAQ();
  StateAQ & operator=(const StateAQ &);

/// Interpolate full fields
  void changeResolution(const StateAQ & xx);

/// Interactions with Increment
  StateAQ & operator+=(const IncrementAQ &);

/// I/O and diagnostics
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const {return fields_->norm();}
  const util::DateTime & validTime() const {return fields_->time();}
  util::DateTime & validTime() {return fields_->time();}
  void updateTime(const util::Duration & dt) {fields_->updateTime(dt);}

/// Access to fields
  FieldsAQ & fields() {return *fields_;}
  const FieldsAQ & fields() const {return *fields_;}
  std::shared_ptr<const GeometryAQ> geometry() const {
    return fields_->geometry();
  }
  const oops::Variables & variables() const {return fields_->variables();}

/// Serialization
  size_t serialSize() const;
  void serialize(std::vector<double> &) const;
  void deserialize(const std::vector<double> &, size_t &);

/// Other
  void zero();
  void accumul(const double &, const StateAQ &);

 private:
  void print(std::ostream &) const;
  std::unique_ptr<FieldsAQ> fields_;
};
// -----------------------------------------------------------------------------

}  // namespace aq

#endif  // AQ_STATEAQ_H_
