/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2017-2019 UCAR.
 * (C) Copyright 2021-2022 CERFACS.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef AQ_INTERFACE_H_
#define AQ_INTERFACE_H_

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace util {
  class DateTime;
  class Duration;
}

namespace aq {
  class Locations;
  class ObsSpace;

// Geometry key type
typedef int F90geom;
// Geometry iterator key type
typedef int F90iter;
// Model key type
typedef int F90model;
// Geovals key type
typedef int F90geovals;
// Fields key type
typedef int F90flds;
// Observation vector key type
typedef int F90ovec;
// Observation data base key type
typedef int F90odb;
// Interpolator key type
typedef int F90interp;

}  // namespace aq
#endif  // AQ_INTERFACE_H_
