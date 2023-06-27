/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2021-2022 CERFACS.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef AQ_MODELDATA_H_
#define AQ_MODELDATA_H_

#include <iostream>
#include <string>

#include "eckit/memory/NonCopyable.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace eckit {
  class LocalConfiguration;
}

namespace aq {
  class Geometry;

/// Model data for the AQ model.
/*!
 * This class is used by SABER/VADER to obtain model-specific constants
 * VADER is not used yet, so this is currently empty
 */

// -----------------------------------------------------------------------------

class ModelData : public util::Printable,
                  private eckit::NonCopyable,
                  private util::ObjectCounter<ModelData> {
 public:
  static const std::string classname() {return "aq::ModelData";}

  explicit ModelData(const Geometry &) {}
  virtual ~ModelData() {}
  ModelData(const ModelData &) = delete;
  ModelData(ModelData &&) = default;
  const ModelData & operator=(const ModelData &) = delete;
  ModelData & operator=(ModelData &&) = default;

  const eckit::LocalConfiguration modelData() const {
    eckit::LocalConfiguration modelData;
    return modelData;
  }

 private:
  void print(std::ostream &) const {}
};

// -----------------------------------------------------------------------------

}  // namespace aq

#endif  // AQ_MODELDATA_H_
