/*
 * (C) Crown Copyright 2021, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "aq/FinalCheck.h"

namespace aq {
// -----------------------------------------------------------------------------
static oops::interface::FilterMaker<AqObsTraits, FinalCheck> makerPreChk_("Final Check");
// -----------------------------------------------------------------------------
}  // namespace aq
