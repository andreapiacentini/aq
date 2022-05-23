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
