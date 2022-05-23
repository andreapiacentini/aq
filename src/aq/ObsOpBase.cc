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

#include "aq/ObsOpBase.h"

#include "eckit/config/Configuration.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

namespace aq {

// -----------------------------------------------------------------------------

ObsOpFactory::ObsOpFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    oops::Log::error() << name << " already registered in aq::ObsOpFactory." << std::endl;
    ABORT("Element already registered in ObsOpFactory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

ObsOpBase * ObsOpFactory::create(const ObsSpace & odb,
                                   const eckit::Configuration & conf) {
  oops::Log::trace() << "ObsOpBase::create starting" << std::endl;
  oops::Log::debug() << "ObsOpBase::create conf" << conf << std::endl;
  const std::string id = conf.getString("obs type");
  typename std::map<std::string, ObsOpFactory*>::iterator
    jloc = getMakers().find(id);
  if (jloc == getMakers().end()) {
    oops::Log::error() << id << " does not exist in observation operator factory." << std::endl;
    ABORT("Element does not exist in ObsOpFactory.");
  }
  ObsOpBase * ptr = jloc->second->make(odb, conf);
  oops::Log::trace() << "ObsOpBase::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

}  // namespace aq
