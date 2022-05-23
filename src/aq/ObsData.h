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

#ifndef AQ_OBSDATA_H_
#define AQ_OBSDATA_H_

#include <math.h>
#include <ostream>
#include <string>
#include <vector>

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "aq/ObsSpace.h"
#include "aq/ObsVec.h"

namespace aq {

// -----------------------------------------------------------------------------
/// Data in observation space

template<typename DATATYPE>
class ObsData : public util::Printable,
                  private util::ObjectCounter<ObsData<DATATYPE> > {
 public:
  static const std::string classname() {return "aq::ObsData";}

  ObsData(const ObsSpace &, const oops::Variables &, const std::string &);
  ObsData(const ObsData &);
  explicit ObsData(const ObsVec &);
  ~ObsData() {}

  ObsData & operator= (const ObsData &);

  /// set all values to zero
  void zero();
  /// set all values to one
  void ones();
  void mask(const ObsData<int>);

// I/O
  void read(const std::string &);
  void save(const std::string &) const;

  const int & toFortran() const {return data_.toFortran();}
  const ObsVec & vect() const {return data_;}
 private:
  void print(std::ostream &) const;

  ObsVec data_;
};
//-----------------------------------------------------------------------------

template<typename DATATYPE>
ObsData<DATATYPE>::ObsData(const ObsSpace & os, const oops::Variables & var,
                               const std::string & name): data_(os, name) {
}
// -----------------------------------------------------------------------------
template<typename DATATYPE>
ObsData<DATATYPE>::ObsData(const ObsData & other): data_(other.data_) {
}
// -----------------------------------------------------------------------------
template<typename DATATYPE>
ObsData<DATATYPE>::ObsData(const ObsVec & other): data_(other) {
}
// -----------------------------------------------------------------------------
template<typename DATATYPE>
ObsData<DATATYPE> & ObsData<DATATYPE>::operator= (const ObsData & rhs) {
  data_ = rhs.data_;
  return *this;
}
// -----------------------------------------------------------------------------
template<typename DATATYPE>
void ObsData<DATATYPE>::zero() {
  data_.zero();
}
// -----------------------------------------------------------------------------
template<typename DATATYPE>
void ObsData<DATATYPE>::ones() {
  data_.ones();
}
// -----------------------------------------------------------------------------
template<typename DATATYPE>
void ObsData<DATATYPE>::mask(const ObsData<int> mask) {
  data_.mask(mask);
}
// -----------------------------------------------------------------------------
template<typename DATATYPE>
void ObsData<DATATYPE>::read(const std::string & name) {
  data_.read(name);
}
// -----------------------------------------------------------------------------
template<typename DATATYPE>
void ObsData<DATATYPE>::save(const std::string & name) const {
  data_.save(name);
}
// -----------------------------------------------------------------------------
template<typename DATATYPE>
void ObsData<DATATYPE>::print(std::ostream & os) const {
  os << data_;
}
// -----------------------------------------------------------------------------
}  // namespace aq

#endif  // AQ_OBSDATA_H_
