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

#ifndef AQ_GEOVALS_H_
#define AQ_GEOVALS_H_

#include <ostream>
#include <string>
#include <vector>

#include "eckit/mpi/Comm.h"

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/Printable.h"

#include "aq/interface.h"

namespace aq {
  class Locations;

/// Parameters controlling a GeoVaLs read/write
class GeoValsParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(GeoValsParameters, Parameters)

 public:
  oops::RequiredParameter<std::string> filename{"filename", "filename for input and output",
                                                this};
};

/// GeoVals class to handle local model values for AQ model.

class GeoVals : public util::Printable,
                private util::ObjectCounter<GeoVals> {
 public:
  typedef GeoValsParameters Parameters_;

  static const std::string classname() {return "aq::GeoVals";}

  GeoVals(const Locations &, const oops::Variables &, const std::vector<size_t> &);
  GeoVals(const Parameters_ &, const ObsSpace &, const oops::Variables &);
  explicit GeoVals(const GeoVals &);

  GeoVals(): keyGeoVals_(0), comm_(eckit::mpi::comm(NULL)) {}
  explicit GeoVals(int & fgeovals): keyGeoVals_(fgeovals), comm_(eckit::mpi::self()) {}

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
  void read(const Parameters_ &);
  void write(const Parameters_ &) const;

  /// communicator
  const eckit::mpi::Comm & comm() const {return comm_;}

  const int & toFortran() const {return keyGeoVals_;}

  void fill(const std::vector<size_t> &, const std::vector<double> &);
  void fillAD(const std::vector<size_t> &, std::vector<double> &) const;

 private:
  void print(std::ostream &) const;
  F90geovals keyGeoVals_;
  oops::Variables vars_;
  const eckit::mpi::Comm & comm_;
};

}  // namespace aq

#endif  // AQ_GEOVALS_H_
