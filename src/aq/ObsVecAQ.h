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

#ifndef AQ_MODEL_OBSVECAQ_H_
#define AQ_MODEL_OBSVECAQ_H_

#include <Eigen/Dense>
#include <ostream>
#include <string>

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "oops/aq/AqFortran.h"

namespace aq {
  class ObsSpaceAQ;
  template <typename DATATYPE> class ObsDataAQ;

// -----------------------------------------------------------------------------
/// ObsVecAQ class to handle vectors in observation space for AQ model.

class ObsVecAQ : public util::Printable,
                 private util::ObjectCounter<ObsVecAQ> {
 public:
  static const std::string classname() {return "aq::ObsVecAQ";}

  ObsVecAQ(const ObsSpaceAQ &,
           const std::string & name = "");
  ObsVecAQ(const ObsVecAQ &);
  ~ObsVecAQ();

  ObsVecAQ & operator = (const ObsVecAQ &);
  ObsVecAQ & operator*= (const double &);
  ObsVecAQ & operator+= (const ObsVecAQ &);
  ObsVecAQ & operator-= (const ObsVecAQ &);
  ObsVecAQ & operator*= (const ObsVecAQ &);
  ObsVecAQ & operator/= (const ObsVecAQ &);

  Eigen::VectorXd packEigen(const ObsVecAQ &) const;
  size_t packEigenSize(const ObsVecAQ &) const;
  size_t size() const;

  /// set all values to zero
  void zero();
  /// set \p i-th value to missing value
  void setToMissing(int i);
  /// set all values to one
  void ones();
  void axpy(const double &, const ObsVecAQ &);
  void invert();
  void random();
  double dot_product_with(const ObsVecAQ &) const;
  double rms() const;
  void mask(const ObsDataAQ<int> &);
  void mask(const ObsVecAQ &);
  ObsVecAQ & operator=(const ObsDataAQ<float> &);

  unsigned int nobs() const;

  const int & toFortran() const {return keyOvec_;}

// I/O
  void save(const std::string &) const;
  void read(const std::string &);

 private:
  void print(std::ostream &) const;

  const ObsSpaceAQ & obsdb_;
  F90ovec keyOvec_;
};
// -----------------------------------------------------------------------------

}  // namespace aq

#endif  // AQ_MODEL_OBSVECAQ_H_
