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

#ifndef AQ_OBSVEC_H_
#define AQ_OBSVEC_H_

#include <Eigen/Dense>
#include <ostream>
#include <string>

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "aq/interface.h"

namespace aq {
  class ObsSpace;
  template <typename DATATYPE> class ObsData;

// -----------------------------------------------------------------------------
/// ObsVec class to handle vectors in observation space for AQ model.

class ObsVec : public util::Printable,
                 private util::ObjectCounter<ObsVec> {
 public:
  static const std::string classname() {return "aq::ObsVec";}

  ObsVec(const ObsSpace &,
         const std::string & name = "", const bool useObservedVariables = false);
  ObsVec(const ObsVec &);
  ~ObsVec();

  ObsVec & operator = (const ObsVec &);
  ObsVec & operator*= (const double &);
  ObsVec & operator+= (const ObsVec &);
  ObsVec & operator-= (const ObsVec &);
  ObsVec & operator*= (const ObsVec &);
  ObsVec & operator/= (const ObsVec &);

  Eigen::VectorXd packEigen(const ObsVec &) const;
  size_t packEigenSize(const ObsVec &) const;
  size_t size() const;

  /// set all values to zero
  void zero();
  /// set \p i-th value to missing value
  void setToMissing(int i);
  /// set all values to one
  void ones();
  void axpy(const double &, const ObsVec &);
  void invert();
  void random();
  double dot_product_with(const ObsVec &) const;
  double rms() const;
  void mask(const ObsData<int> &);
  void mask(const ObsVec &);
  ObsVec & operator=(const ObsData<float> &);

  unsigned int nobs() const;

  const int & toFortran() const {return keyOvec_;}

// I/O
  void save(const std::string &) const;
  void read(const std::string &);

 private:
  void print(std::ostream &) const;
  const eckit::mpi::Comm & comm_;
  const ObsSpace & obsdb_;
  F90ovec keyOvec_;
};
// -----------------------------------------------------------------------------

}  // namespace aq

#endif  // AQ_OBSVEC_H_
