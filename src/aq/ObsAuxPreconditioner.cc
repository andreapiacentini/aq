/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <vector>

#include "aq/ObsAuxPreconditioner.h"

namespace aq {
  class ObsAuxControl;
// -----------------------------------------------------------------------------

ObsAuxPreconditioner::ObsAuxPreconditioner(const std::array<double,
                                           ObsAuxControl::ntypes> & precond) :
    precond_(precond) {}

void ObsAuxPreconditioner::multiply(const ObsAuxIncrement & dx1, ObsAuxIncrement & dx2) const {
    for (unsigned int jj = 0; jj < ObsAuxControl::ntypes; ++jj) {
      if (precond_[jj] > 0.0) {
        dx2[jj] = dx1[jj] * precond_[jj];
      } else {
        dx2[jj] = 0.0;
      }
    }
}

// -----------------------------------------------------------------------------

}  // namespace aq
