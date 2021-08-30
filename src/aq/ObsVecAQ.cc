/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2017-2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <math.h>

#include "oops/util/Logger.h"

#include "model/AqFortran.h"
#include "model/ObsDataAQ.h"
#include "model/ObsSpaceAQ.h"
#include "model/ObsVecAQ.h"

#include "eckit/exception/Exceptions.h"

namespace aq {
// -----------------------------------------------------------------------------
ObsVecAQ::ObsVecAQ(const ObsSpaceAQ & obsdb, const std::string & name)
  : obsdb_(obsdb), keyOvec_(0)
{
  aq_obsvec_setup_f90(keyOvec_, obsdb.obsvariables().size(), obsdb.nobs());
  if (!name.empty()) obsdb_.getdb(name, keyOvec_);
}
// -----------------------------------------------------------------------------
ObsVecAQ::ObsVecAQ(const ObsVecAQ & other)
  : obsdb_(other.obsdb_), keyOvec_(0) {
  aq_obsvec_clone_f90(keyOvec_, other.keyOvec_);
  aq_obsvec_copy_f90(keyOvec_, other.keyOvec_);
}
// -----------------------------------------------------------------------------
ObsVecAQ::~ObsVecAQ() {
  aq_obsvec_delete_f90(keyOvec_);
}
// -----------------------------------------------------------------------------
ObsVecAQ & ObsVecAQ::operator= (const ObsVecAQ & rhs) {
  const int keyOvecRhs = rhs.keyOvec_;
  aq_obsvec_copy_f90(keyOvec_, keyOvecRhs);
  return *this;
}
// -----------------------------------------------------------------------------
ObsVecAQ & ObsVecAQ::operator*= (const double & zz) {
  aq_obsvec_mul_scal_f90(keyOvec_, zz);
  return *this;
}
// -----------------------------------------------------------------------------
ObsVecAQ & ObsVecAQ::operator+= (const ObsVecAQ & rhs) {
  const int keyOvecRhs = rhs.keyOvec_;
  aq_obsvec_add_f90(keyOvec_, keyOvecRhs);
  return *this;
}
// -----------------------------------------------------------------------------
ObsVecAQ & ObsVecAQ::operator-= (const ObsVecAQ & rhs) {
  const int keyOvecRhs = rhs.keyOvec_;
  aq_obsvec_sub_f90(keyOvec_, keyOvecRhs);
  return *this;
}
// -----------------------------------------------------------------------------
ObsVecAQ & ObsVecAQ::operator*= (const ObsVecAQ & rhs) {
  const int keyOvecRhs = rhs.keyOvec_;
  aq_obsvec_mul_f90(keyOvec_, keyOvecRhs);
  return *this;
}
// -----------------------------------------------------------------------------
ObsVecAQ & ObsVecAQ::operator/= (const ObsVecAQ & rhs) {
  const int keyOvecRhs = rhs.keyOvec_;
  aq_obsvec_div_f90(keyOvec_, keyOvecRhs);
  return *this;
}
// -----------------------------------------------------------------------------
ObsVecAQ & ObsVecAQ::operator=(const ObsDataAQ<float> & rhs) {
  *this = rhs.vect();
  return *this;
}
// -----------------------------------------------------------------------------
void ObsVecAQ::zero() {
  aq_obsvec_zero_f90(keyOvec_);
}
// -----------------------------------------------------------------------------
void ObsVecAQ::setToMissing(int ii) {
  aq_obsvec_settomissing_ith_f90(keyOvec_, ii);
}
// -----------------------------------------------------------------------------
void ObsVecAQ::ones() {
  aq_obsvec_ones_f90(keyOvec_);
}
// -----------------------------------------------------------------------------
void ObsVecAQ::axpy(const double & zz, const ObsVecAQ & rhs) {
  const int keyOvecRhs = rhs.keyOvec_;
  aq_obsvec_axpy_f90(keyOvec_, zz, keyOvecRhs);
}
// -----------------------------------------------------------------------------
void ObsVecAQ::invert() {
  aq_obsvec_invert_f90(keyOvec_);
}
// -----------------------------------------------------------------------------
void ObsVecAQ::random() {
  aq_obsvec_random_f90(obsdb_, keyOvec_);
}
// -----------------------------------------------------------------------------
double ObsVecAQ::dot_product_with(const ObsVecAQ & other) const {
  const int keyOvecOther = other.keyOvec_;
  double zz;
  aq_obsvec_dotprod_f90(keyOvec_, keyOvecOther, zz);
  return zz;
}
// -----------------------------------------------------------------------------
double ObsVecAQ::rms() const {
  int iobs;
  aq_obsvec_nobs_f90(keyOvec_, iobs);
  double zz = 0.0;
  if (iobs > 0) {
    aq_obsvec_dotprod_f90(keyOvec_, keyOvec_, zz);
    zz = sqrt(zz/iobs);
  }
  return zz;
}
// -----------------------------------------------------------------------------
void ObsVecAQ::mask(const ObsDataAQ<int> & mask) {
  aq_obsvec_mask_f90(keyOvec_, mask.toFortran());
}
// -----------------------------------------------------------------------------
void ObsVecAQ::mask(const ObsVecAQ & mask) {
  aq_obsvec_mask_with_missing_f90(keyOvec_, mask.toFortran());
}
// -----------------------------------------------------------------------------
void ObsVecAQ::save(const std::string & name) const {
  obsdb_.putdb(name, keyOvec_);
}
// -----------------------------------------------------------------------------
Eigen::VectorXd ObsVecAQ::packEigen(const ObsVecAQ & mask) const {
  Eigen::VectorXd vec(packEigenSize(mask));
  aq_obsvec_get_withmask_f90(keyOvec_, mask.toFortran(), vec.data(), vec.size());
  return vec;
}
// -----------------------------------------------------------------------------
size_t ObsVecAQ::packEigenSize(const ObsVecAQ & mask) const {
  int nobs;
  aq_obsvec_nobs_withmask_f90(keyOvec_, mask.toFortran(), nobs);
  return nobs;
}
// -----------------------------------------------------------------------------
void ObsVecAQ::read(const std::string & name) {
  obsdb_.getdb(name, keyOvec_);
}
// -----------------------------------------------------------------------------
void ObsVecAQ::print(std::ostream & os) const {
  if (nobs() == 0) {
    os << obsdb_.obsname() << " no observations.";
  } else {
    double zmin, zmax, zavg;
    aq_obsvec_stats_f90(keyOvec_, zmin, zmax, zavg);
    std::ios_base::fmtflags f(os.flags());
    os << obsdb_.obsname() << " nobs= " << nobs()
       << std::scientific << std::setprecision(4)
       << "  Min=" << std::setw(12) << zmin
       << ", Max=" << std::setw(12) << zmax
       << ", Average=" << std::setw(12) << zavg;
    os.flags(f);
  }
}
// -----------------------------------------------------------------------------
unsigned int ObsVecAQ::nobs() const {
  int iobs;
  aq_obsvec_nobs_f90(keyOvec_, iobs);
  unsigned int nobs(iobs);
  return nobs;
}
// -----------------------------------------------------------------------------
size_t ObsVecAQ::size() const {
  int iobs;
  aq_obsvec_size_f90(keyOvec_, iobs);
  size_t nobs(iobs);
  return nobs;
}
// -----------------------------------------------------------------------------
}  // namespace aq
