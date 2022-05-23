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

#include <mpi.h>

#include <fstream>

#include "eckit/utils/Translator.h"

#include "aq/Tools.h"

#include "oops/mpi/mpi.h"

namespace aq {
  // ---------------------------------------------------------------------
  // getEnv utility from LibOOPS
  //
  bool getEnv(const std::string& env, bool default_value) {
    if (::getenv(env.c_str()))
      {return eckit::Translator<std::string, bool>()(::getenv(env.c_str()));}
    return default_value;
  }

  int getEnv(const std::string& env, int default_value) {
    if (::getenv(env.c_str()))
      {return eckit::Translator<std::string, int>()(::getenv(env.c_str()));}
    return default_value;
  }

  void syncAll() {
    oops::mpi::world().barrier();
  }

  std::ofstream ofnull_;

  std::ostream& masterOut() {
    if (eckit::mpi::comm().rank() == 0) {
      return std::cout;
    } else {
      ofnull_.setstate(std::ios_base::badbit);
      return ofnull_;
    }
  }

  int timeInstance() {
    return static_cast<int>(oops::mpi::world().rank() / eckit::mpi::comm().size());
  }

  int nbTimeInstances() {
    return static_cast<int>(oops::mpi::world().size() / eckit::mpi::comm().size());
  }

  // -----------------------------------------------------------------------------

}  // namespace aq
