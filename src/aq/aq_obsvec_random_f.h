/*
 * (C) Copyright 2017-2020 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef AQ_MODEL_AQ_OBSVEC_RANDOM_F_H_
#define AQ_MODEL_AQ_OBSVEC_RANDOM_F_H_

namespace aq {
  class ObsSpaceAQ;

extern "C" {
  void aq_obsvec_random_f(const ObsSpaceAQ &, const int &, double *);
}

}  // namespace aq

#endif  // AQ_MODEL_AQ_OBSVEC_RANDOM_F_H_
