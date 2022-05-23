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
