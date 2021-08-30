/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef AQ_FIELDSAQ_H_
#define AQ_FIELDSAQ_H_

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

#include "aq/GeometryAQ.h"
#include "aq/GeometryAQIterator.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
  class LocalIncrement;
}

namespace aq {
  class LocationsAQ;
  class GomAQ;

// -----------------------------------------------------------------------------
/// Class to represent a Fields for the AQ model
class FieldsAQ : public util::Printable,
                 public util::Serializable,
                 private util::ObjectCounter<FieldsAQ> {
 public:
  static const std::string classname() {return "aq::FieldsAQ";}

// Constructors and basic operators
  FieldsAQ(const GeometryAQ &, const oops::Variables &, const util::DateTime &);
  FieldsAQ(const FieldsAQ &, const GeometryAQ &);  // Not implemented
  FieldsAQ(const FieldsAQ &, const oops::Variables &);
  FieldsAQ(const FieldsAQ &, const bool);
  FieldsAQ(const FieldsAQ &);
  ~FieldsAQ();

  void zero();
  void zero(const util::DateTime &);
  void ones();
  FieldsAQ & operator=(const FieldsAQ &);
  FieldsAQ & operator+=(const FieldsAQ &);
  FieldsAQ & operator-=(const FieldsAQ &);
  FieldsAQ & operator*=(const double &);
  void axpy(const double &, const FieldsAQ &);
  double dot_product_with(const FieldsAQ &) const;
  void schur_product_with(const FieldsAQ &);
  void dirac(const eckit::Configuration &);
  void random();

// Interpolate full fields
  void changeResolution(const FieldsAQ &);
  void add(const FieldsAQ &);
  void diff(const FieldsAQ &, const FieldsAQ &);

// ATLAS FieldSet
  void setAtlas(atlas::FieldSet *) const;
  void toAtlas(atlas::FieldSet *) const;
  void fromAtlas(atlas::FieldSet *);

// Utilities
  void read(const eckit::Configuration &);
  void analytic_init(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const;
  std::shared_ptr<const GeometryAQ> geometry() const {return geom_;}
  const oops::Variables & variables() const {return vars_;}

  const util::DateTime & time() const {return time_;}
  util::DateTime & time() {return time_;}
  void updateTime(const util::Duration & dt) {time_ += dt;}

  const int & toFortran() const {return keyFlds_;}

  bool isForModel(const bool &) const {}

  oops::LocalIncrement getLocal(const GeometryAQIterator &) const;
  void setLocal(const oops::LocalIncrement &, const GeometryAQIterator &);

/// Serialization
  size_t serialSize() const override;
  void serialize(std::vector<double> &) const override;
  void deserialize(const std::vector<double> &, size_t &) override;

 private:
  void print(std::ostream &) const override;
  F90flds keyFlds_;
  std::shared_ptr<const GeometryAQ> geom_;
  const oops::Variables vars_;
  util::DateTime time_;
};
// -----------------------------------------------------------------------------

}  // namespace aq
#endif  // AQ_FIELDSAQ_H_
