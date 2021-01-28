/* -------------------------------------------------------------------------- *
 *                     OpenSim:  HaeufleActiveForceLengthCurve.cpp            *
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

#include "HaeufleActiveForceLengthCurve.h"

using namespace OpenSim;
using namespace SimTK;
using namespace std;

//==============================================================================
// CONSTRUCTION
//==============================================================================
// Uses default (compiler-generated) destructor, copy constructor, and copy
// assignment.

// Default constructor.
HaeufleActiveForceLengthCurve::HaeufleActiveForceLengthCurve() {
    setNull();
    constructProperties();
    setName(getConcreteClassName());
}

HaeufleActiveForceLengthCurve::HaeufleActiveForceLengthCurve(
    double exponentDescendingActiveForceLength,
    double widthDescendingActiveForceLength,
    double exponentAscendingActiveForceLength,
    double widthAscendingActiveForceLength) {
    
    setNull();
    constructProperties();
    setName(getConcreteClassName());

    set_exponent_descending_active_force_length(
            exponentDescendingActiveForceLength);
    set_width_ascending_active_force_length(widthDescendingActiveForceLength);
    set_exponent_ascending_active_force_length(
            exponentAscendingActiveForceLength);
    set_width_ascending_active_force_length(widthAscendingActiveForceLength);
}

void HaeufleActiveForceLengthCurve::setNull() { setAuthors("Mike Spahr"); }

void HaeufleActiveForceLengthCurve::constructProperties() 
{ 
    constructProperty_exponent_descending_active_force_length(1.50);
    constructProperty_width_ascending_active_force_length(0.45);
    constructProperty_exponent_ascending_active_force_length(3.00);
    constructProperty_width_ascending_active_force_length(0.45);

}

//==============================================================================
// COMPONENT INTERFACE
//==============================================================================
void HaeufleActiveForceLengthCurve::extendFinalizeFromProperties() {
    Super::extendFinalizeFromProperties();

    std::string errorLocation =
            getName() +
            " HaeufleActiveForceLengthCurve::extendFinalizeFromProperties";

    // Ensure property values are within appropriate ranges.
    OPENSIM_THROW_IF_FRMOBJ(get_exponent_descending_active_force_length() <= 0,
            InvalidPropertyValue,
            getProperty_exponent_descending_active_force_length().getName(),
            "The exponent of the descending branch of the active force length "
            "must be greater than zero");
    OPENSIM_THROW_IF_FRMOBJ(get_width_descending_active_force_length() <= 0,
            InvalidPropertyValue,
            getProperty_width_descending_active_force_length().getName(),
            "The width of the descending branch of the active force length "
            "curve must be greater than zero.");
    OPENSIM_THROW_IF_FRMOBJ(get_exponent_ascending_active_force_length() <= 0,
            InvalidPropertyValue,
            getProperty_exponent_ascending_active_force_length().getName(),
            "The exponent of the ascending branch of the active force length "
            "must be greater than zero.");
    OPENSIM_THROW_IF_FRMOBJ(get_width_ascending_active_force_length() <= 0,
            InvalidPropertyValue,
            getProperty_width_ascending_active_force_length().getName(),
            "The width of the ascending branch of the active force length must "
            "be greater than zero.");
}


/** @returns The exponent of the descending limb of the normalized
    bell curve */
double
HaeufleActiveForceLengthCurve::getExponentDescendingActiveForceLength() const {
    return get_exponent_descending_active_force_length();
}

/** @returns The width of the descending limb of the normalized
bell curve */
double
HaeufleActiveForceLengthCurve::getWidthDescendingActiveForceLength() const {
    return get_width_descending_active_force_length();
}

/** @returns The exponent of the ascending limb of the normalized
bell curve */
double
HaeufleActiveForceLengthCurve::getExponentAscendingActiveForceLength() const {
    return get_exponent_ascending_active_force_length();
}

/** @returns The width of the ascending limb of the normalized
bell curve */
double
HaeufleActiveForceLengthCurve::getWidthAscendingActiveForceLength() const {
    return get_width_ascending_active_force_length();
}

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
void HaeufleActiveForceLengthCurve::setActiveFiberLengths(
    double exponentDescendingActiveForceLength,
    double widthDescendingActiveForceLength,
    double exponentAscendingActiveForceLength,
    double widthAscendingActiveForceLength) {
    set_exponent_descending_active_force_length(
            exponentDescendingActiveForceLength);
    set_width_descending_active_force_length(widthDescendingActiveForceLength);
    set_exponent_ascending_active_force_length(
            exponentAscendingActiveForceLength);
    set_width_ascending_active_force_length(widthAscendingActiveForceLength);
}

/** Evaluates the active-force-length curve at a normalized fiber length of
'normFiberLength'. */
double HaeufleActiveForceLengthCurve::calcValue(double normFiberLength) const {
    // define variables and set them to values which will let the simulation fail 
    // if they are not overwritten.
    double exponent_active_force_length = 0;
    double width_active_force_length = 0;
    // if normalized fiber length is greater than 1 we are in the descending domain
    if (normFiberLength > 1) {
        exponent_active_force_length = getExponentDescendingActiveForceLength();
        width_active_force_length = getWidthDescendingActiveForceLength();
    } else { // if normalized fiber length is smaller than 1 we are in the ascending domain
        exponent_active_force_length = getExponentAscendingActiveForceLength();
        width_active_force_length = getWidthAscendingActiveForceLength();
    }
    return exp(-pow(abs((normFiberLength - 1) / width_active_force_length),
                    exponent_active_force_length));
}
