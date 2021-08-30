/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2017-2019 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef AQ_INCREMENTAQ_H_
#define AQ_INCREMENTAQ_H_

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

#include "aq/FieldsAQ.h"
#include "aq/GeometryAQ.h"
#include "aq/GeometryAQIterator.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class LocalIncrement;
  class Variables;
}

namespace aq {
  class GomAQ;
  class LocationsAQ;
  class GeometryAQ;
  class ModelBiasIncrement;
  class ErrorCovarianceAQ;
  class StateAQ;

/// Increment Class: Difference between two states
/*!
 *  Some fields that are present in a State may not be present in
 *  an Increment. The Increment contains everything that is needed by
 *  the tangent-linear and adjoint models.
 */

// -----------------------------------------------------------------------------

class IncrementAQ : public util::Printable,
                    public util::Serializable,
                    private util::ObjectCounter<IncrementAQ> {
 public:
  static const std::string classname() {return "aq::IncrementAQ";}

/// Constructor, destructor
  IncrementAQ(const GeometryAQ &, const oops::Variables &, const util::DateTime &);
  IncrementAQ(const GeometryAQ &, const IncrementAQ &);
  IncrementAQ(const IncrementAQ &, const bool);
  IncrementAQ(const IncrementAQ &);
  virtual ~IncrementAQ();

/// Basic operators
  void diff(const StateAQ &, const StateAQ &);
  void zero();
  void zero(const util::DateTime &);
  void ones();
  IncrementAQ & operator =(const IncrementAQ &);
  IncrementAQ & operator+=(const IncrementAQ &);
  IncrementAQ & operator-=(const IncrementAQ &);
  IncrementAQ & operator*=(const double &);
  void axpy(const double &, const IncrementAQ &, const bool check = true);
  double dot_product_with(const IncrementAQ &) const;
  void schur_product_with(const IncrementAQ &);
  void random();
  void dirac(const eckit::Configuration &);

/// I/O and diagnostics
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const {return fields_->norm();}
  const util::DateTime & validTime() const {return fields_->time();}
  util::DateTime & validTime() {return fields_->time();}
  void updateTime(const util::Duration & dt) {fields_->time() += dt;}

/// ATLAS FieldSet
  void setAtlas(atlas::FieldSet *) const;
  void toAtlas(atlas::FieldSet *) const;
  void fromAtlas(atlas::FieldSet *);

/// Access to fields
  FieldsAQ & fields() {return *fields_;}
  const FieldsAQ & fields() const {return *fields_;}

  std::shared_ptr<const GeometryAQ> geometry() const {
    return fields_->geometry();
  }

/// Other
  void accumul(const double &, const StateAQ &);
  oops::LocalIncrement getLocal(const GeometryAQIterator &) const;
  void setLocal(const oops::LocalIncrement &, const GeometryAQIterator &);

/// Serialization
  size_t serialSize() const override;
  void serialize(std::vector<double> &) const override;
  void deserialize(const std::vector<double> &, size_t &) override;

/// Data
 private:
  void print(std::ostream &) const override;
  std::unique_ptr<FieldsAQ> fields_;
};
// -----------------------------------------------------------------------------

}  // namespace aq

#endif  // AQ_INCREMENTAQ_H_
