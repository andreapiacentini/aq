/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef AQ_GOMAQ_H_
#define AQ_GOMAQ_H_

#include <ostream>
#include <string>
#include <vector>

#include "oops/base/Variables.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "aq/AqFortran.h"

namespace oops {
  class Variables;
}

namespace aq {
  class LocationsAQ;

/// GomAQ class to handle local model values for AQ model.

class GomAQ : public util::Printable,
              private util::ObjectCounter<GomAQ> {
 public:
  static const std::string classname() {return "aq::GomAQ";}

  GomAQ(const LocationsAQ &, const oops::Variables &, const std::vector<size_t> &);
  GomAQ(const eckit::Configuration &, const ObsSpaceAQ &,
        const oops::Variables &);
  explicit GomAQ(const GomAQ &);

  GomAQ(): keyGom_(0) {}
  explicit GomAQ(int & fgom): keyGom_(fgom) {}

  ~GomAQ();

  void zero();
  void random();
  double rms() const;
  double normalizedrms(const GomAQ &) const;
  GomAQ & operator=(const GomAQ &);
  GomAQ & operator*=(const double &);
  GomAQ & operator+=(const GomAQ &);
  GomAQ & operator-=(const GomAQ &);
  GomAQ & operator*=(const GomAQ &);
  double dot_product_with(const GomAQ &) const;
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;

  const int & toFortran() const {return keyGom_;}

 private:
  void print(std::ostream &) const;
  F90gom keyGom_;
  oops::Variables vars_;
};

}  // namespace aq

#endif  // AQ_GOMAQ_H_
