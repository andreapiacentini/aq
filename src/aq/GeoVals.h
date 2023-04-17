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

#ifndef AQ_GEOVALS_H_
#define AQ_GEOVALS_H_

#include <Eigen/Core>

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
  /// References to read-only or writable vector- or matrix-valued expressions.
  ///
  /// For example, an Eigen::Vector, Eigen::Matrix or an Eigen::Map (the latter can be used as a
  /// view onto a chunk of memory stored in another container, such as a std::vector).
  template <typename T>
  using ConstVectorRef = Eigen::Ref<const Eigen::Vector<T, Eigen::Dynamic>>;
  template <typename T>
  using ConstMatrixRef = Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>;
  template <typename T>
  using MatrixRef = Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>;

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

  void fill(const std::string &name, const ConstVectorRef<size_t> &indx,
            const ConstMatrixRef<double> &vals, const bool levelsTopDown);
  void fillAD(const std::string &name, const ConstVectorRef<size_t> &indx,
              MatrixRef<double> vals, const bool levelsTopDown) const;

 private:
  void print(std::ostream &) const;
  F90geovals keyGeoVals_;
  oops::Variables vars_;
  const eckit::mpi::Comm & comm_;
};

}  // namespace aq

#endif  // AQ_GEOVALS_H_
