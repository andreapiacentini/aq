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

#ifndef AQ_FIELDS_H_
#define AQ_FIELDS_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/LocalIncrement.h"
#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Serializable.h"

#include "aq/Geometry.h"
#include "aq/GeometryIterator.h"
#include "aq/interface.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
  class LocalIncrement;
}

namespace aq {
  class Locations;
  class GeoVals;

// -----------------------------------------------------------------------------
/// Class to represent a Fields for the AQ model
class Fields : public util::Printable,
                 public util::Serializable,
                 private util::ObjectCounter<Fields> {
 public:
  static const std::string classname() {return "aq::Fields";}

// Constructors and basic operators
  Fields(const Geometry &, const oops::Variables &, const util::DateTime &);
  Fields(const Fields &, const Geometry &);  // Not implemented
  Fields(const Fields &, const oops::Variables &);
  Fields(const Fields &, const bool);
  Fields(const Fields &);
  ~Fields();

  void zero();
  void zero(const util::DateTime &);
  void ones();
  Fields & operator=(const Fields &);
  Fields & operator+=(const Fields &);
  Fields & operator-=(const Fields &);
  Fields & operator*=(const double &);
  void axpy(const double &, const Fields &);
  double dot_product_with(const Fields &) const;
  void schur_product_with(const Fields &);
  void dirac(const eckit::Configuration &);
  void random();

// Interpolate full fields
  void changeResolution(const Fields &);
  void add(const Fields &);
  void diff(const Fields &, const Fields &);

// ATLAS FieldSet
  void setAtlas(atlas::FieldSet *) const;
  void toAtlas(atlas::FieldSet *) const;
  void fromAtlas(atlas::FieldSet *);

// Utilities
  void read(const eckit::Configuration &);
  void analytic_init(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const;
  std::shared_ptr<const Geometry> geometry() const {return geom_;}
  const oops::Variables & variables() const {return vars_;}

  const util::DateTime & time() const {return time_;}
  util::DateTime & time() {return time_;}
  void updateTime(const util::Duration & dt) {time_ += dt;}

  const int & toFortran() const {return keyFlds_;}

  oops::LocalIncrement getLocal(const GeometryIterator &) const;
  void setLocal(const oops::LocalIncrement &, const GeometryIterator &);

/// Serialization
  size_t serialSize() const override;
  void serialize(std::vector<double> &) const override;
  void deserialize(const std::vector<double> &, size_t &) override;

 private:
  void print(std::ostream &) const override;
  F90flds keyFlds_;
  std::shared_ptr<const Geometry> geom_;
  const oops::Variables vars_;
  util::DateTime time_;
};
// -----------------------------------------------------------------------------

}  // namespace aq
#endif  // AQ_FIELDS_H_
