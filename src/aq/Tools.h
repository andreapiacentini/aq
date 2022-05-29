/*
 * (C) Copyright 2019-2020 UCAR
 * (C) Copyright 2021-2022 CERFACS.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef AQ_TOOLS_H_
#define AQ_TOOLS_H_

#include <string>

namespace aq {

  // ---------------------------------------------------------------------
  // getEnv utility from LibOOPS
  //
  bool getEnv(const std::string& env, bool default_value);

  int getEnv(const std::string& env, int default_value);

  void syncAll();

  std::ostream& masterOut();

  int timeInstance();

  int nbTimeInstances();

  // -----------------------------------------------------------------------------

}  // namespace aq

#endif  // AQ_TOOLS_H_
