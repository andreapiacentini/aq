/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef AQ_MODELAUXINCREMENT_H_
#define AQ_MODELAUXINCREMENT_H_

#include <iostream>
#include <vector>

#include "oops/util/Printable.h"
#include "oops/util/Serializable.h"

namespace eckit {
  class Configuration;
}

namespace aq {
  class ModelAuxControl;
  class ModelAuxCovariance;
  class Geometry;

// -----------------------------------------------------------------------------

class ModelAuxIncrement : public util::Printable,
                           public util::Serializable {
 public:
/// Constructor, destructor
  ModelAuxIncrement(const Geometry &, const eckit::Configuration &) {}
  ModelAuxIncrement(const ModelAuxIncrement &, const bool) {}
  ModelAuxIncrement(const ModelAuxIncrement &, const eckit::Configuration &) {}
  ~ModelAuxIncrement() {}

/// Linear algebra operators
  void diff(const ModelAuxControl &, const ModelAuxControl &) {}
  void zero() {}
  ModelAuxIncrement & operator=(const ModelAuxIncrement &) {return *this;}
  ModelAuxIncrement & operator+=(const ModelAuxIncrement &) {return *this;}
  ModelAuxIncrement & operator-=(const ModelAuxIncrement &) {return *this;}
  ModelAuxIncrement & operator*=(const double) {return *this;}
  void axpy(const double, const ModelAuxIncrement &) {}
  double dot_product_with(const ModelAuxIncrement &) const {return 0.0;}

/// I/O and diagnostics
  void read(const eckit::Configuration &) {}
  void write(const eckit::Configuration &) const {}
  double norm() const {return 0.0;}

/// Serialization
  size_t serialSize() const override {return 0;}
  void serialize(std::vector<double> &) const override {}
  void deserialize(const std::vector<double> &, size_t &) override {}

 private:
  explicit ModelAuxIncrement(const ModelAuxCovariance &);
  void print(std::ostream & os) const override {}
};

// -----------------------------------------------------------------------------

}  // namespace aq

#endif  // AQ_MODELAUXINCREMENT_H_
