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

#ifndef AQ_INCREMENT_H_
#define AQ_INCREMENT_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/LocalIncrement.h"
#include "oops/util/DateTime.h"
#include "oops/util/dot_product.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Serializable.h"

#include "aq/Fields.h"
#include "aq/Geometry.h"
#include "aq/GeometryIterator.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class LocalIncrement;
  class Variables;
}

namespace aq {
  class GeoVals;
  class Locations;
  class Geometry;
  class ModelAuxIncrement;
  class Covariance;
  class State;

/// Increment Class: Difference between two states
/*!
 *  Some fields that are present in a State may not be present in
 *  an Increment. The Increment contains everything that is needed by
 *  the tangent-linear and adjoint models.
 */

// -----------------------------------------------------------------------------

class Increment : public util::Printable,
                    public util::Serializable,
                    private util::ObjectCounter<Increment> {
 public:
  static const std::string classname() {return "aq::Increment";}

/// Constructor, destructor
  Increment(const Geometry &, const oops::Variables &, const util::DateTime &);
  Increment(const Geometry &, const Increment &);
  Increment(const Increment &, const bool);
  Increment(const Increment &);
  virtual ~Increment();

/// Basic operators
  void diff(const State &, const State &);
  void zero();
  void zero(const util::DateTime &);
  void ones();
  Increment & operator =(const Increment &);
  Increment & operator+=(const Increment &);
  Increment & operator-=(const Increment &);
  Increment & operator*=(const double &);
  void axpy(const double &, const Increment &, const bool check = true);
  double dot_product_with(const Increment &) const;
  void schur_product_with(const Increment &);
  void random();
  void dirac(const eckit::Configuration &);

/// I/O and diagnostics
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const {return fields_->norm();}
  std::vector<double> rmsByLevel(const std::string &) const {}
  const util::DateTime & validTime() const {return fields_->time();}
  util::DateTime & validTime() {return fields_->time();}
  void updateTime(const util::Duration & dt) {fields_->time() += dt;}

/// ATLAS FieldSet
  void toFieldSet(atlas::FieldSet &) const;
  void toFieldSetAD(const atlas::FieldSet &);
  void fromFieldSet(const atlas::FieldSet &);

/// Access to fields
  Fields & fields() {return *fields_;}
  const Fields & fields() const {return *fields_;}

/// Other
  void accumul(const double &, const State &);
  oops::LocalIncrement getLocal(const GeometryIterator &) const;
  void setLocal(const oops::LocalIncrement &, const GeometryIterator &);

/// Serialization
  size_t serialSize() const override;
  void serialize(std::vector<double> &) const override;
  void deserialize(const std::vector<double> &, size_t &) override;

/// Variables
  const oops::Variables & variables() const {return vars_;}

/// Data
 private:
  void print(std::ostream &) const override;
  std::unique_ptr<Fields> fields_;
  oops::Variables vars_;
};
// -----------------------------------------------------------------------------

}  // namespace aq

#endif  // AQ_INCREMENT_H_
