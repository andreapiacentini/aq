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
