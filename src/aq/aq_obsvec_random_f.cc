/*
 * (C) Copyright 2017-2020 UCAR
 * (C) Copyright 2021-2022 CERFACS.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "aq/aq_obsvec_random_f.h"
#include "aq/ObsSpace.h"
#include "oops/util/Random.h"

// -----------------------------------------------------------------------------
namespace aq {
// -----------------------------------------------------------------------------

void aq_obsvec_random_f(const ObsSpace & odb, const int & nn, double * xx) {
  // AQ  static util::NormalDistribution<double> dist(nn, 0.0, 1.0, odb.getSeed());
  util::NormalDistribution<double> dist(nn, 0.0, 1.0, odb.getSeed());
  for (int jj = 0; jj < nn; ++jj) xx[jj] = dist[jj];
}

// -----------------------------------------------------------------------------

}  // namespace aq
