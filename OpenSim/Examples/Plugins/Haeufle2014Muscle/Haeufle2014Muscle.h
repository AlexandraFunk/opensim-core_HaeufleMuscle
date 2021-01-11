#ifndef OPENSIM_Haeufle2014Muscle_h__
#define OPENSIM_Haeufle2014Muscle_h__
/* -------------------------------------------------------------------------- *
 *                  OpenSim:  Haeufle2014Muscle.h                  *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Copyright (c) 2005-2017 Stanford University and the Authors                *
 * Author(s): Matthew Millard, Tom Uchida, Ajay Seth                          *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */
#include "osimPluginDLL.h"
#include <simbody/internal/common.h>

// The parent class, Muscle.h, provides
//    1. max_isometric_force
//    2. optimal_fiber_length
//    3. tendon_slack_length
//    4. pennation_angle_at_optimal
//    5. max_contraction_velocity
//    6. ignore_tendon_compliance
//    7. ignore_activation_dynamics
#include <OpenSim/Simulation/Model/Muscle.h>

// Sub-models used by this muscle model
// TODO add Rockenfeller activation dynamic model
#include <OpenSim/Actuators/MuscleFirstOrderActivationDynamicModel.h>
#include <OpenSim/Actuators/MuscleFixedWidthPennationModel.h>

#include "HaeufleActiveForceLengthCurve.h"
#include "HaeufleForceVelocityCurve.h"
#include "HaeufleFiberForceLengthCurve.h"
#include "HaeufleTendonForceLengthCurve.h"
#include "HaeufleTendonDampingCurve.h"
// #include <OpenSim/Actuators/HaeufleParallelDampingCurve.h> // never used nor modelled in demoa before 

#ifdef SWIG
    #ifdef OSIMACTUATORS_API
        #undef OSIMACTUATORS_API
        #define OSIMACTUATORS_API
    #endif
#endif

namespace OpenSim {

//==============================================================================
//                         Haeufle2014Muscle
//==============================================================================
/**
This class implements a configurable muscle model, as described in
Haufle et al.\ (2013). 

TODO: Add more description

<B>Reference</B>

D.F.B. Haeufle, M. Guenther, A. Bayer, S. Schmitt (2014) Hill-type muscle model 
with serial damping and eccentric force–velocity relation. Journal of
Biomechanics https://doi.org/10.1016/j.jbiomech.2014.02.009.

@author Maria Hammer
@author Mike Spahr
*/

} //end of namespace OpenSim

#endif // OPENSIM_Haeufle2014Muscle_h__
