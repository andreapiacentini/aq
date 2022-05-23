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

#include <math.h>

#include "oops/util/Logger.h"

#include "aq/aq_obsvec_interface.h"
#include "aq/ObsData.h"
#include "aq/ObsSpace.h"
#include "aq/ObsVec.h"

#include "eckit/exception/Exceptions.h"

namespace aq {
// -----------------------------------------------------------------------------
  ObsVec::ObsVec(const ObsSpace & obsdb, const std::string & name,
                 const bool useObservedVariables)
  : obsdb_(obsdb), keyOvec_(0), comm_(obsdb.comm())
  {
    if (comm_.rank() == 0) {
      aq_obsvec_setup_f90(keyOvec_, obsdb.assimvariables().size(), obsdb.nobs());
      if (!name.empty()) obsdb_.getdb(name, keyOvec_);
    }
  }
// -----------------------------------------------------------------------------
  ObsVec::ObsVec(const ObsVec & other)
  : obsdb_(other.obsdb_), keyOvec_(0), comm_(other.comm_) {
    if (comm_.rank() == 0) {
      aq_obsvec_clone_f90(keyOvec_, other.keyOvec_);
      aq_obsvec_copy_f90(keyOvec_, other.keyOvec_);
    }
  }
// -----------------------------------------------------------------------------
  ObsVec::~ObsVec() {
    if (comm_.rank() == 0) {
      aq_obsvec_delete_f90(keyOvec_);
    }
  }
// -----------------------------------------------------------------------------
  ObsVec & ObsVec::operator= (const ObsVec & rhs) {
    if (comm_.rank() == 0) {
      const int keyOvecRhs = rhs.keyOvec_;
      aq_obsvec_copy_f90(keyOvec_, keyOvecRhs);
    }
    return *this;
  }
// -----------------------------------------------------------------------------
  ObsVec & ObsVec::operator*= (const double & zz) {
    if (comm_.rank() == 0) {
      aq_obsvec_mul_scal_f90(keyOvec_, zz);
    }
    return *this;
  }
// -----------------------------------------------------------------------------
  ObsVec & ObsVec::operator+= (const ObsVec & rhs) {
    if (comm_.rank() == 0) {
      const int keyOvecRhs = rhs.keyOvec_;
      aq_obsvec_add_f90(keyOvec_, keyOvecRhs);
    }
    return *this;
  }
// -----------------------------------------------------------------------------
  ObsVec & ObsVec::operator-= (const ObsVec & rhs) {
    if (comm_.rank() == 0) {
      const int keyOvecRhs = rhs.keyOvec_;
      aq_obsvec_sub_f90(keyOvec_, keyOvecRhs);
    }
    return *this;
  }
// -----------------------------------------------------------------------------
  ObsVec & ObsVec::operator*= (const ObsVec & rhs) {
    if (comm_.rank() == 0) {
      const int keyOvecRhs = rhs.keyOvec_;
      aq_obsvec_mul_f90(keyOvec_, keyOvecRhs);
    }
    return *this;
  }
// -----------------------------------------------------------------------------
  ObsVec & ObsVec::operator/= (const ObsVec & rhs) {
    if (comm_.rank() == 0) {
      const int keyOvecRhs = rhs.keyOvec_;
      aq_obsvec_div_f90(keyOvec_, keyOvecRhs);
    }
    return *this;
  }
// -----------------------------------------------------------------------------
  ObsVec & ObsVec::operator=(const ObsData<float> & rhs) {
    if (comm_.rank() == 0) {
      *this = rhs.vect();
    }
    return *this;
  }
// -----------------------------------------------------------------------------
  void ObsVec::zero() {
    if (comm_.rank() == 0) {
      aq_obsvec_zero_f90(keyOvec_);
    }
  }
// -----------------------------------------------------------------------------
  void ObsVec::setToMissing(int ii) {
    if (comm_.rank() == 0) {
      aq_obsvec_settomissing_ith_f90(keyOvec_, ii);
    }
  }
// -----------------------------------------------------------------------------
  void ObsVec::ones() {
    if (comm_.rank() == 0) {
      aq_obsvec_ones_f90(keyOvec_);
    }
  }
// -----------------------------------------------------------------------------
  void ObsVec::axpy(const double & zz, const ObsVec & rhs) {
    if (comm_.rank() == 0) {
      const int keyOvecRhs = rhs.keyOvec_;
      aq_obsvec_axpy_f90(keyOvec_, zz, keyOvecRhs);
    }
  }
// -----------------------------------------------------------------------------
  void ObsVec::invert() {
    if (comm_.rank() == 0) {
      aq_obsvec_invert_f90(keyOvec_);
    }
  }
// -----------------------------------------------------------------------------
  void ObsVec::random() {
    if (comm_.rank() == 0) {
      aq_obsvec_random_f90(obsdb_, keyOvec_);
    }
  }
// -----------------------------------------------------------------------------
  double ObsVec::dot_product_with(const ObsVec & other) const {
    double zz;
    if (comm_.rank() == 0) {
      const int keyOvecOther = other.keyOvec_;
      aq_obsvec_dotprod_f90(keyOvec_, keyOvecOther, zz);
    }
    comm_.broadcast(zz, 0);
    return zz;
  }
// -----------------------------------------------------------------------------
  double ObsVec::rms() const {
    double zz = 0.0;
    if (comm_.rank() == 0) {
      int iobs;
      aq_obsvec_nobs_f90(keyOvec_, iobs);
      if (iobs > 0) {
        aq_obsvec_dotprod_f90(keyOvec_, keyOvec_, zz);
        zz = sqrt(zz/iobs);
      }
    }
    comm_.broadcast(zz, 0);
    return zz;
  }
// -----------------------------------------------------------------------------
  void ObsVec::mask(const ObsData<int> & mask) {
    if (comm_.rank() == 0) {
      aq_obsvec_mask_f90(keyOvec_, mask.toFortran());
    }
  }
// -----------------------------------------------------------------------------
  void ObsVec::mask(const ObsVec & mask) {
    if (comm_.rank() == 0) {
      aq_obsvec_mask_with_missing_f90(keyOvec_, mask.toFortran());
    }
  }
// -----------------------------------------------------------------------------
  void ObsVec::save(const std::string & name) const {
    if (comm_.rank() == 0) {
      obsdb_.putdb(name, keyOvec_);
    }
  }
// -----------------------------------------------------------------------------
  Eigen::VectorXd ObsVec::packEigen(const ObsVec & mask) const {
    Eigen::VectorXd vec(packEigenSize(mask));
    if (comm_.rank() == 0) {
      aq_obsvec_get_withmask_f90(keyOvec_, mask.toFortran(), vec.data(), vec.size());
    }
    return vec;
  }
// -----------------------------------------------------------------------------
  size_t ObsVec::packEigenSize(const ObsVec & mask) const {
    int nobs;
    if (comm_.rank() == 0) {
      aq_obsvec_nobs_withmask_f90(keyOvec_, mask.toFortran(), nobs);
    }
    comm_.broadcast(nobs, 0);
    return nobs;
  }
// -----------------------------------------------------------------------------
  void ObsVec::read(const std::string & name) {
    if (comm_.rank() == 0) {
      obsdb_.getdb(name, keyOvec_);
    }
  }
// -----------------------------------------------------------------------------
  void ObsVec::print(std::ostream & os) const {
    unsigned int iobs = nobs();
    if (comm_.rank() == 0) {
      if (iobs == 0) {
        os << obsdb_.obsname() << " no observations.";
      } else {
        double zmin, zmax, zavg;
        aq_obsvec_stats_f90(keyOvec_, zmin, zmax, zavg);
        std::ios_base::fmtflags f(os.flags());
        os << obsdb_.obsname() << " nobs= " << iobs
        << std::scientific << std::setprecision(4)
        << "  Min=" << std::setw(12) << zmin
        << ", Max=" << std::setw(12) << zmax
        << ", Average=" << std::setw(12) << zavg;
        os.flags(f);
      }
    }
  }
// -----------------------------------------------------------------------------
  unsigned int ObsVec::nobs() const {
    int iobs;
    if (comm_.rank() == 0) {
      aq_obsvec_nobs_f90(keyOvec_, iobs);
    }
    comm_.broadcast(iobs, 0);
    unsigned int nobs(iobs);
    return nobs;
  }
// -----------------------------------------------------------------------------
  size_t ObsVec::size() const {
    int iobs;
    if (comm_.rank() == 0) {
      aq_obsvec_size_f90(keyOvec_, iobs);
    }
    comm_.broadcast(iobs, 0);
    size_t nobs(iobs);
    return nobs;
  }
// -----------------------------------------------------------------------------
}  // namespace aq
