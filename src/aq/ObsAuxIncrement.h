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
