/* -------------------------------------------------------------------------- *
 *                    OpenSim:  HaeufleTendonForceLengthCurve.cpp             *
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
#include "HaeufleTendonForceLengthCurve.h"

using namespace OpenSim;
using namespace SimTK;
using namespace std;

HaeufleTendonForceLengthCurve::HaeufleTendonForceLengthCurve() 
{
    setNull();
    constructProperties();
    setName(getConcreteClassName());
}

HaeufleTendonForceLengthCurve::HaeufleTendonForceLengthCurve(
        double tendonSlackLength, double forceAtNonlinearLinearTransition)
{
    setNull();
    constructProperties();
    setName(getConcreteClassName());

    set_tendon_slack_length(tendonSlackLength);
    set_force_at_nonlinear_linear_transition(forceAtNonlinearLinearTransition);

}

HaeufleTendonForceLengthCurve::HaeufleTendonForceLengthCurve(
        double tendonSlackLength,
    double relativeStretchAtNonlinearLinearTransition,
    double relativeStretchAtLinearPart,
    double forceAtNonlinearLinearTransition) 
{
    setNull();
    constructProperties();
    setName(getConcreteClassName());

    set_tendon_slack_length(tendonSlackLength);
    set_relative_stretch_at_linear_part(relativeStretchAtLinearPart);
    set_relative_stretch_at_nonlinear_linear_transition(
            relativeStretchAtNonlinearLinearTransition);
    set_force_at_nonlinear_linear_transition(forceAtNonlinearLinearTransition);
}

void HaeufleTendonForceLengthCurve::setNull() 
{ setAuthors("Mike Spahr"); }

void HaeufleTendonForceLengthCurve::constructProperties() 
{
    constructProperty_tendon_slack_length(0.15);
    constructProperty_relative_stretch_at_nonlinear_linear_transition(0.0425);
    constructProperty_relative_stretch_at_linear_part(0.0170);
    constructProperty_force_at_nonlinear_linear_transition(500);
}



void HaeufleTendonForceLengthCurve::extendFinalizeFromProperties() {
    Super::extendFinalizeFromProperties();

    std::string errorLocation =
            getName() +
            " HaeufleTendonForceLengthCurve::extendFinalizeFromProperties";

    // Ensure property values are within appropriate ranges.
    OPENSIM_THROW_IF_FRMOBJ(get_tendon_slack_length() <= 0,
            InvalidPropertyValue, getProperty_tendon_slack_length().getName(),
            "The serial elastic rest length must be greater than zero");
    OPENSIM_THROW_IF_FRMOBJ(
            get_relative_stretch_at_nonlinear_linear_transition() <= 0,
            InvalidPropertyValue,
            getProperty_relative_stretch_at_nonlinear_linear_transition()
                    .getName(),
            "The relative stretch at nonlinear/linear transition must be "
            "greater than zero");
    OPENSIM_THROW_IF_FRMOBJ(get_relative_stretch_at_linear_part() <= 0,
            InvalidPropertyValue,
            getProperty_relative_stretch_at_linear_part().getName(),
            "The relative stretch in the linear part must be greater than "
            "zero");
    OPENSIM_THROW_IF_FRMOBJ(get_force_at_nonlinear_linear_transition() <= 0,
            InvalidPropertyValue,
            getProperty_force_at_nonlinear_linear_transition().getName(),
            "The force at the nonlinear/linear transition must be greater than "
            "zero");
}

double HaeufleTendonForceLengthCurve::getTendonSlackLength() const {
    return get_tendon_slack_length();
}

double
HaeufleTendonForceLengthCurve::getRelativeStretchAtNonlinearLinearTransition()
const {
    return get_relative_stretch_at_nonlinear_linear_transition();

}


double HaeufleTendonForceLengthCurve::getRelativeStretchAtLinearPart() const {
    return get_relative_stretch_at_linear_part();
}

double
HaeufleTendonForceLengthCurve::getForceAtNonlinearLinearTransition() const {
    return get_force_at_nonlinear_linear_transition();
}

void HaeufleTendonForceLengthCurve::setTendonSlackLength(
        double aTendonSlackLength) {
    set_tendon_slack_length(aTendonSlackLength);
}


void HaeufleTendonForceLengthCurve::setForceAtNonlinearLinearTransition(
    double aForceAtNonlinearLinearTransition) {
    set_force_at_nonlinear_linear_transition(aForceAtNonlinearLinearTransition);
}

void HaeufleTendonForceLengthCurve::setRelativeStretchAtNonlinearLinearTransition(
    double aRelativeStretchAtNonlinearLinearTransition) {
    set_relative_stretch_at_nonlinear_linear_transition(
            aRelativeStretchAtNonlinearLinearTransition);
}


void HaeufleTendonForceLengthCurve::setRelativeStretchAtLinearPart(
    double aRelativeStretchAtLinearPart) {
    set_relative_stretch_at_linear_part(aRelativeStretchAtLinearPart);
}

double HaeufleTendonForceLengthCurve::calcValue(
        double aSerialElasticLength) const {
    double Lsee0 = getTendonSlackLength();
    double Fsee =
            0; // initialize Fsee for first case aSerialElasticLength < Lsee0
    double deltaUseenll = getRelativeStretchAtNonlinearLinearTransition();
    if (aSerialElasticLength > Lsee0 &&
            aSerialElasticLength <= (Lsee0 * (1 + deltaUseenll))) 
    {
        double deltaFsee0 = getForceAtNonlinearLinearTransition();
        double nuesee = deltaUseenll / getRelativeStretchAtLinearPart();
        Fsee = deltaFsee0 *
               pow((aSerialElasticLength - Lsee0) / (deltaUseenll * Lsee0),
                       nuesee);
    } 
    else if (aSerialElasticLength > (Lsee0*(1+deltaUseenll))) 
    {
        double deltaFsee0 = getForceAtNonlinearLinearTransition();
        double deltaUseel = getRelativeStretchAtLinearPart();
        Fsee = deltaFsee0 *
               (1 + (aSerialElasticLength - Lsee0 * (1 + deltaUseenll)) /
                               (deltaUseel * Lsee0));
    }
    return Fsee;
}

