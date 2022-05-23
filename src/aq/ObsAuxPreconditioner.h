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

#ifndef AQ_OBSAUXPRECONDITIONER_H_
#define AQ_OBSAUXPRECONDITIONER_H_

#include <ostream>
#include <string>
#include <boost/noncopyable.hpp>

#include "aq/ObsAuxControl.h"
#include "aq/ObsAuxIncrement.h"
#include "aq/ObsAuxParameters.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace aq {
// -----------------------------------------------------------------------------

/// Class for VarBC preconditioner.

class ObsAuxPreconditioner : public util::Printable,
                             private boost::noncopyable,
                             private util::ObjectCounter<ObsAuxPreconditioner> {
 public:
  typedef ObsAuxControlParameters Parameters_;

  static const std::string classname() {return "aq::ObsAuxPreconditioner";}

/// Constructor, destructor
  explicit ObsAuxPreconditioner(const std::array<double, ObsAuxControl::ntypes>  &);

/// Linear algebra operators
  void multiply(const ObsAuxIncrement & dx1, ObsAuxIncrement & dx2) const;

 private:
  void print(std::ostream &) const {}
  const std::array<double, ObsAuxControl::ntypes> precond_;
};

// -----------------------------------------------------------------------------

}  // namespace aq

#endif  // AQ_OBSAUXPRECONDITIONER_H_
