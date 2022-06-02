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
  const Geometry & geom_;
  const oops::Variables vars_;
  util::DateTime time_;
};
// -----------------------------------------------------------------------------

}  // namespace aq
#endif  // AQ_FIELDS_H_
