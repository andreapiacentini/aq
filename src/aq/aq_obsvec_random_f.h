/*
 * (C) Copyright 2017-2020 UCAR
 * (C) Copyright 2021-2022 CERFACS.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef AQ_AQ_OBSVEC_RANDOM_F_H_
#define AQ_AQ_OBSVEC_RANDOM_F_H_

namespace aq {
  class ObsSpace;

extern "C" {
  void aq_obsvec_random_f(const ObsSpace &, const int &, double *);
}

}  // namespace aq

#endif  // AQ_AQ_OBSVEC_RANDOM_F_H_
