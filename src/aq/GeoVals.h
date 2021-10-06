/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef AQ_GEOVALS_H_
#define AQ_GEOVALS_H_

#include <ostream>
#include <string>
#include <vector>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "aq/interface.h"

namespace oops {
  class Variables;
}

namespace aq {
  class Locations;

/// GeoVals class to handle local model values for AQ model.

class GeoVals : public util::Printable,
              private util::ObjectCounter<GeoVals> {
 public:
  static const std::string classname() {return "aq::GeoVals";}

  GeoVals(const Locations &, const oops::Variables &, const std::vector<size_t> &);
  GeoVals(const eckit::Configuration &, const ObsSpace &,
        const oops::Variables &);
  explicit GeoVals(const GeoVals &);

  GeoVals(): keyGeoVals_(0) {}
  explicit GeoVals(int & fgeovals): keyGeoVals_(fgeovals) {}

  ~GeoVals();

  void zero();
  void random();
  double rms() const;
  double normalizedrms(const GeoVals &) const;
  GeoVals & operator=(const GeoVals &);
  GeoVals & operator*=(const double &);
  GeoVals & operator+=(const GeoVals &);
  GeoVals & operator-=(const GeoVals &);
  GeoVals & operator*=(const GeoVals &);
  double dot_product_with(const GeoVals &) const;
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;

  const int & toFortran() const {return keyGeoVals_;}

 private:
  void print(std::ostream &) const;
  F90geovals keyGeoVals_;
  oops::Variables vars_;
};

}  // namespace aq

#endif  // AQ_GEOVALS_H_