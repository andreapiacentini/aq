/*
 * (C) Copyright 2018-2019 UCAR
 * (C) Copyright 2021-2022 CERFACS.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "aq/FinalCheck.h"

namespace aq {
// -----------------------------------------------------------------------------
static oops::interface::FilterMaker<ObsTraits, FinalCheck> makerQCm_("Final Check");
// -----------------------------------------------------------------------------
}  // namespace aq
