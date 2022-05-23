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

#ifndef AQ_OBSAUXINCREMENT_H_
#define AQ_OBSAUXINCREMENT_H_

#include <iostream>
#include <vector>

#include "oops/util/Printable.h"
#include "oops/util/Serializable.h"

#include "aq/ObsAuxParameters.h"

namespace eckit {
  class Configuration;
}

namespace aq {
  class ObsAuxControl;
  class ObsSpace;

// -----------------------------------------------------------------------------

class ObsAuxIncrement : public util::Printable,
                         public util::Serializable {
 public:
  typedef ObsAuxControlParameters Parameters_;

/// Constructor, destructor
  ObsAuxIncrement();
  ObsAuxIncrement(const ObsSpace &, const Parameters_ &);
  ObsAuxIncrement(const ObsAuxIncrement &, const bool copy = true);

/// Linear algebra operators
  void diff(const ObsAuxControl &, const ObsAuxControl &);
  void zero();
  ObsAuxIncrement & operator=(const ObsAuxIncrement &);
  ObsAuxIncrement & operator+=(const ObsAuxIncrement &);
  ObsAuxIncrement & operator-=(const ObsAuxIncrement &);
  ObsAuxIncrement & operator*=(const double);
  void axpy(const double, const ObsAuxIncrement &);
  double dot_product_with(const ObsAuxIncrement &) const;

/// I/O and diagnostics
  void read(const eckit::Configuration &) {}
  void write(const eckit::Configuration &) const {}
  double norm() const;

  double & operator[](const unsigned int ii) {return bias_[ii];}
  const double & operator[](const unsigned int ii) const {return bias_[ii];}

  double & insitu() {return bias_[0];}
  const double & insitu() const {return bias_[0];}

/// Serialization
  size_t serialSize() const override;
  void serialize(std::vector<double> &) const override;
  void deserialize(const std::vector<double> &, size_t &) override;

 private:
  void print(std::ostream &) const override;
  void makePassive();

  std::vector<double> bias_;
  std::vector<bool> active_;
};

// -----------------------------------------------------------------------------

}  // namespace aq

#endif  // AQ_OBSAUXINCREMENT_H_
