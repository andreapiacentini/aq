/*
 * (C) Copyright 2017-2018 UCAR
 * (C) Copyright 2021-2022 CERFACS.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "aq/LinearObsOpBase.h"

#include "eckit/config/Configuration.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

namespace aq {

// -----------------------------------------------------------------------------

LinearObsOpFactory::LinearObsOpFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    oops::Log::error() << name << " already registered in aq::LinearObsOpFactory." << std::endl;
    ABORT("Element already registered in LinearObsOpFactory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

LinearObsOpBase * LinearObsOpFactory::create(const ObsSpace & odb,
                                         const eckit::Configuration & conf) {
  oops::Log::trace() << "ObsOpBase::create starting" << std::endl;
  oops::Log::debug() << "ObsOpBase::create conf" << conf << std::endl;
  const std::string id = conf.getString("obs type");
  typename std::map<std::string, LinearObsOpFactory*>::iterator
    jloc = getMakers().find(id);
  if (jloc == getMakers().end()) {
    oops::Log::error() << id << " does not exist in observation operator factory." << std::endl;
    ABORT("Element does not exist in LinearObsOpFactory.");
  }
  LinearObsOpBase * ptr = jloc->second->make(odb, conf);
  oops::Log::trace() << "ObsOpBase::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

}  // namespace aq
