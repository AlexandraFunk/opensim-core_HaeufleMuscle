#ifndef OPENSIM_HAEUFLE_ACTIVE_FORCE_LENGTH_CURVE_H_
#define OPENSIM_HAEUFLE_ACTIVE_FORCE_LENGTH_CURVE_H_


/* -------------------------------------------------------------------------- *
 *                     OpenSim:  HaeufleActiveForceLengthCurve.h              *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Copyright (c) 2005-2017 Stanford University and the Authors                *
 * Author(s): Mike Spahr                                                      *
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

// INCLUDE
#include <OpenSim/Actuators/osimActuatorsDLL.h>
#include <OpenSim/Simulation/Model/ModelComponent.h>

#ifdef SWIG
#	ifdef OSIMACTUATORS_API
#		undef OSIMACTUATORS_API
#		define OSIMACTUATORS_API
#	endif
#endif

namespace OpenSim {
/** This class serves as a ActiveForceLengthCurve as described by Haeufle et 
al in Hill-type muscle model with serial damping and eccentric force–velocity 
relation. Journal of Biomechanics (D.F.B. Haeufle, M. Guenther, A. Bayer, S. 
Schmitt (2014)) 
In contrast to the quintic bezier splines, which is used in the 
ActiveForceLenghtCurve, this implementation uses an exponential function with
different parameters for the ascending and descending part of the graph.
The provided fiber length must be normalized with the optimalFiberLength

@param exponentDescendingActiveForceLength
    The exponent of the normalized bell curve in its descending limb as used in 
    Haeufle et al (2014).
@param widthDescendingActiveForceLength
    The width of the normalized bell curve in its descending limb as used in 
    Haeufle et al (2014).
@param exponentAscendingActiveForceLength
    The exponent of the normalized bell curve in its ascending limb as used in
    Haeufle et al (2014).
@param widthAscendingActiveForceLength
    The width of the normalized bell curve in its ascending limb as used in
    Haeufle et al (2014).

<B>Default Parameter Values</B>
The default parameter values have been chosen to match the files which are 
included in the scripts/dataligament files of calcman. All parameters are 
in SI-Units.
\verbatim
exponentDescendingActiveForceLength .... 1.50
widthDescendingActiveForceLength ....... 0.45
exponentAscendingActiveForceLength ..... 3.00
widthAscendingActiveForceLength ........ 0.45
\endverbatim

<B>Example</B>
@code
HaeufleActiveForceLengthCurve fisomCurve1(1.5, 0.45, 3.00, 0.45);
double fisomVal  = fisomCurve1.calcValue(1.0);
@endcode

<B>References</B>
\li D.F.B. Haeufle, M. Guenther, A. Bayer, S. Schmitt(2014) Hill-type 
    muscle model with serial damping and eccentric force–velocity 
    relation. Journal of Biomechanics

@author Mike Spahr
*/
class OSIMACTUATORS_API HaeufleActiveForceLengthCurve : public ModelComponent {
    OpenSim_DECLARE_CONCRETE_OBJECT(
            HaeufleActiveForceLengthCurve, ModelComponent);

public:
    //==============================================================================
    // PROPERTIES
    //==============================================================================
    OpenSim_DECLARE_PROPERTY(exponent_descending_active_force_length, double,
            "The exponent of the normalized bell curve in its descending"
            "limb as used in Haeufle et al (2014).");
    OpenSim_DECLARE_PROPERTY(width_descending_active_force_length, double,
            "The width of the normalized bell curve in its descending"
            "limb as used in Haeufle et al (2014).");
    OpenSim_DECLARE_PROPERTY(exponent_ascending_active_force_length, double,
            "The exponent of the normalized bell curve in its ascending"
            "limb as used in Haeufle et al (2014).");
    OpenSim_DECLARE_PROPERTY(width_ascending_active_force_length, double,
            "The width of the normalized bell curve in its ascending"
            "limb as used in Haeufle et al (2014).");

    //==============================================================================
    // PUBLIC METHODS
    //==============================================================================
    /** The default constructor creates an active-force-length curve using the
    default property values and assigns a default name. */
    HaeufleActiveForceLengthCurve();

    /** Constructs an active-force-length curve using the provided parameters
    and assigns a default name. */
    HaeufleActiveForceLengthCurve(double exponentDescendingActiveForceLength,
            double widthDescendingActiveForceLength,
            double exponentAscendingActiveForceLength,
            double widthAscendingActiveForceLength);
    
    /** @returns The exponent of the descending limb of the normalized
    bell curve */
    double getExponentDescendingActiveForceLength() const;

    /** @returns The width of the descending limb of the normalized
    bell curve */
    double getWidthDescendingActiveForceLength() const;

    /** @returns The exponent of the ascending limb of the normalized
    bell curve */
    double getExponentAscendingActiveForceLength() const;

    /** @returns The width of the ascending limb of the normalized
    bell curve */
    double getWidthAscendingActiveForceLength() const;
    
    /**
    @param exponentDescendingActiveForceLength
        The exponent of the descending limb of the normalized
        bell curve.
    @param widthDescendingActiveForceLength
        The width of the descending limb of the normalized
        bell curve.
    @param exponentAscendingActiveForceLength
        The exponent of the ascending limb of the normalized
        bell curve.
    @param widthAscendingActiveForceLength
        The width of the ascending limb of the normalized
        bell curve.
    */
    void setActiveFiberLengths(double exponentDescendingActiveForceLength,
            double widthDescendingActiveForceLength,
            double exponentAscendingActiveForceLength,
            double widthAscendingActiveForceLength);

    /** Evaluates the active-force-length curve at a normalized fiber length of
    'normFiberLength'. */
    double calcValue(double normFiberLength) const;

    //==============================================================================
    // PRIVATE
    //==============================================================================
private:
    void setNull();
    void constructProperties();
};

} // namespace OpenSim

#endif // OPENSIM_HAEUFLE_ACTIVE_FORCE_LENGTH_CURVE_H__
