/*
 * (C) Copyright 2019  UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef AQ_MODEL_OBSDATAAQ_H_
#define AQ_MODEL_OBSDATAAQ_H_

#include <math.h>
#include <ostream>
#include <string>
#include <vector>

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "oops/aq/AqFortran.h"
#include "oops/aq/ObsSpaceAQ.h"
#include "oops/aq/ObsVecAQ.h"

namespace aq {

// -----------------------------------------------------------------------------
/// Data in observation space

template<typename DATATYPE>
class ObsDataAQ : public util::Printable,
                  private util::ObjectCounter<ObsDataAQ<DATATYPE> > {
 public:
  static const std::string classname() {return "aq::ObsDataAQ";}

  ObsDataAQ(const ObsSpaceAQ &, const oops::Variables &, const std::string &);
  ObsDataAQ(const ObsDataAQ &);
  explicit ObsDataAQ(const ObsVecAQ &);
  ~ObsDataAQ() {}

  ObsDataAQ & operator= (const ObsDataAQ &);

  /// set all values to zero
  void zero();
  /// set all values to one
  void ones();
  void mask(const ObsDataAQ<int>);

// I/O
  void read(const std::string &);
  void save(const std::string &) const;

  const int & toFortran() const {return data_.toFortran();}
  const ObsVecAQ & vect() const {return data_;}
 private:
  void print(std::ostream &) const;

  ObsVecAQ data_;
};
//-----------------------------------------------------------------------------

template<typename DATATYPE>
ObsDataAQ<DATATYPE>::ObsDataAQ(const ObsSpaceAQ & os, const oops::Variables & var,
                               const std::string & name): data_(os, name) {
}
// -----------------------------------------------------------------------------
template<typename DATATYPE>
ObsDataAQ<DATATYPE>::ObsDataAQ(const ObsDataAQ & other): data_(other.data_) {
}
// -----------------------------------------------------------------------------
template<typename DATATYPE>
ObsDataAQ<DATATYPE>::ObsDataAQ(const ObsVecAQ & other): data_(other) {
}
// -----------------------------------------------------------------------------
template<typename DATATYPE>
ObsDataAQ<DATATYPE> & ObsDataAQ<DATATYPE>::operator= (const ObsDataAQ & rhs) {
  data_ = rhs.data_;
  return *this;
}
// -----------------------------------------------------------------------------
template<typename DATATYPE>
void ObsDataAQ<DATATYPE>::zero() {
  data_.zero();
}
// -----------------------------------------------------------------------------
template<typename DATATYPE>
void ObsDataAQ<DATATYPE>::ones() {
  data_.ones();
}
// -----------------------------------------------------------------------------
template<typename DATATYPE>
void ObsDataAQ<DATATYPE>::mask(const ObsDataAQ<int> mask) {
  aq_obsvec_mask_f90(data_.toFortran(), mask.toFortran());
}
// -----------------------------------------------------------------------------
template<typename DATATYPE>
void ObsDataAQ<DATATYPE>::read(const std::string & name) {
  data_.read(name);
}
// -----------------------------------------------------------------------------
template<typename DATATYPE>
void ObsDataAQ<DATATYPE>::save(const std::string & name) const {
  data_.save(name);
}
// -----------------------------------------------------------------------------
template<typename DATATYPE>
void ObsDataAQ<DATATYPE>::print(std::ostream & os) const {
  os << data_;
}
// -----------------------------------------------------------------------------
}  // namespace aq

#endif  // AQ_MODEL_OBSDATAAQ_H_
