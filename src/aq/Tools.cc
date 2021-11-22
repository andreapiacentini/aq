/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
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
    MPI::COMM_WORLD.Barrier();
  }

  std::ofstream ofnull_;

  std::ostream& masterOut() {
    if (oops::mpi::world().rank() == 0) {
      return std::cout;
    } else {
      ofnull_.setstate(std::ios_base::badbit);
      return ofnull_;
    }
  }

  int timeInstance() {
    return static_cast<int>(MPI::COMM_WORLD.Get_rank() / oops::mpi::world().size());
  }

  int nbTimeInstances() {
    return static_cast<int>(MPI::COMM_WORLD.Get_size() / oops::mpi::world().size());
  }

  // -----------------------------------------------------------------------------

}  // namespace aq
