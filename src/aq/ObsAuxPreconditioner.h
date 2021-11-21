/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
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
