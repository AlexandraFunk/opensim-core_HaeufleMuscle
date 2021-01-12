/* -------------------------------------------------------------------------- *
 *                 OpenSim:  Haeufle2014Muscle.cpp                            *
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
#include "Haeufle2014Muscle.h"
#include <OpenSim/Simulation/Model/Model.h>

using namespace std;
using namespace OpenSim;
using namespace SimTK;

const string Haeufle2014Muscle::STATE_CALCIUM_CONCENTRATION =
        "calcium_concentration";
const string Haeufle2014Muscle::STATE_FIBER_LENGTH_NAME = "fiber_length";

Haeufle2014Muscle::Haeufle2014Muscle() 
{ 
    setNull(); 
    constructProperties();
}


Haeufle2014Muscle::Haeufle2014Muscle(const std::string& aName,
    double aMaxIsometricForce,
    double aOptimalFiberLength, double aTendonSlackLength,
    double aPennationAngle)
{

}

void Haeufle2014Muscle::setNull() { setAuthors("Maria Hammer, Mike Spahr"); }
void Haeufle2014Muscle::constructProperties() {
    constructProperty_default_activation(0.05);
    constructProperty_default_fiber_length(getOptimalFiberLength());

    constructProperty_HaeufleActiveForceLengthCurve(
            HaeufleActiveForceLengthCurve());
    constructProperty_HaeufleForceVelocityCurve(HaeufleForceVelocityCurve());
    constructProperty_HaeufleFiberForceLengthCurve(
            HaeufleFiberForceLengthCurve());
    constructProperty_HaeufleTendonForceLengthCurve(
            HaeufleTendonForceLengthCurve());
    constructProperty_HaeufleTendonDampingCurve(HaeufleTendonDampingCurve());

    //TODO check if this is necessary?
    // setMinControl(get_minimum_activation());

}

void Haeufle2014Muscle::extendFinalizeFromProperties() {
    Super::extendFinalizeFromProperties();

    // TODO add variable dependent properties 
    // throw errors if variables are in the wrong range 
    // check millard2012EquilibriumMuscle function for more info
}


void Haeufle2014Muscle::extendConnectToModel(Model& model) {

}

void Haeufle2014Muscle::extendAddToSystem(
    SimTK::MultibodySystem& system) const {

}

void Haeufle2014Muscle::extendInitStateFromProperties(SimTK::State& s) const {

}

void Haeufle2014Muscle::extendSetPropertiesFromState(const SimTK::State& s) {

}

void Haeufle2014Muscle::computeStateVariableDerivatives(
    const SimTK::State& s) const {

}

//==============================================================================
// GET METHODS
//==============================================================================

double Haeufle2014Muscle::getDefaultFiberLength() const {
    return get_default_fiber_length();
}

const HaeufleActiveForceLengthCurve&
Haeufle2014Muscle::getActiveForceLengthCurve() const {
    return get_HaeufleActiveForceLengthCurve();
}

const HaeufleForceVelocityCurve&
Haeufle2014Muscle::getForceVelocityCurve() const {
    return get_HaeufleForceVelocityCurve();
}

const HaeufleFiberForceLengthCurve&
Haeufle2014Muscle::getFiberForceLengthCurve() const {
    return get_HaeufleFiberForceLengthCurve();
}

const HaeufleTendonForceLengthCurve&
Haeufle2014Muscle::getTendonForceLengthCurve() const {
    return get_HaeufleTendonForceLengthCurve();
}

const HaeufleTendonDampingCurve&
Haeufle2014Muscle::getTendonDampingCurve() const {
    return get_HaeufleTendonDampingCurve();
}

const MuscleFixedWidthPennationModel&
Haeufle2014Muscle::getPennationModel() const {
    return getMemberSubcomponent<MuscleFixedWidthPennationModel>(penMdlIdx);
}

const RockenfellerFirstOrderActivationDynamicModel&
Haeufle2014Muscle::getActivationModel() const {
    return getMemberSubcomponent<RockenfellerFirstOrderActivationDynamicModel>(
            actMdlIdx);
}

//==============================================================================
// SET METHODS
//==============================================================================

void Haeufle2014Muscle::setDefaultActivation(double activation) {
    set_default_activation(activation);
}

void Haeufle2014Muscle::setDefaultFiberLength(double fiberLength) {
    set_default_fiber_length(fiberLength);
}

void Haeufle2014Muscle::setActiveForceLengthCurve(
    HaeufleActiveForceLengthCurve& aActiveForceLengthCurve) {
    set_HaeufleActiveForceLengthCurve(aActiveForceLengthCurve);
}

void Haeufle2014Muscle::setForceVelocityCurve(
    HaeufleForceVelocityCurve& aForceVelocityCurve) {
    set_HaeufleForceVelocityCurve(aForceVelocityCurve);
}

void Haeufle2014Muscle::setFiberForceLengthCurve(
    HaeufleFiberForceLengthCurve& aFiberForceLengthCurve) {
    set_HaeufleFiberForceLengthCurve(aFiberForceLengthCurve);
}

void Haeufle2014Muscle::setTendonForceLengthCurve(
    HaeufleTendonForceLengthCurve& aTendonForceLengthCurve) {
    set_HaeufleTendonForceLengthCurve(aTendonForceLengthCurve);
}

void Haeufle2014Muscle::setTendonDampingCurve(
    HaeufleTendonDampingCurve& aTendonDampingCurve) {
    set_HaeufleTendonDampingCurve(aTendonDampingCurve);
}

void Haeufle2014Muscle::setFiberLength(
    SimTK::State& s, double fiberLength) const {
    if (!get_ignore_tendon_compliance()) {
        setStateVariableValue(
                s, STATE_FIBER_LENGTH_NAME, clampFiberLength(fiberLength));
        markCacheVariableInvalid(s, _lengthInfoCV);
        markCacheVariableInvalid(s, _velInfoCV);
        markCacheVariableInvalid(s, _dynamicsInfoCV);
    }
}

//==============================================================================
// MUSCLE.H INTERFACE
//==============================================================================
void Haeufle2014Muscle::setActivation(
    SimTK::State& s, double activation) const {

}

double Haeufle2014Muscle::computeActuation(const SimTK::State& s) const {

}

void Haeufle2014Muscle::computeFiberEquilibrium(
    SimTK::State& s, bool solveForVelocity = false) const {

}

double Haeufle2014Muscle::getStateVariableDeriv(
    const SimTK::State& s, const std::string& aStateName) const {

}

void Haeufle2014Muscle::setStateVariableDeriv(const SimTK::State& s,
    const std::string& aStateName,
    double aValue) const
{

}

const MuscleLengthInfo& Haeufle2014Muscle::getMuscleLengthInfo(
    const SimTK::State& s) const {

}

MuscleLengthInfo& Haeufle2014Muscle::updMuscleLengthInfo(
    const SimTK::State& s) const {

}

const FiberVelocityInfo& Haeufle2014Muscle::getFiberVelocityInfo(
    const SimTK::State& s) const {

}

FiberVelocityInfo& Haeufle2014Muscle::updFiberVelocityInfo(
    const SimTK::State& s) const
{

}

const MuscleDynamicsInfo& Haeufle2014Muscle::getMuscleDynamicsInfo(
    const SimTK::State& s) const {

}

MuscleDynamicsInfo& Haeufle2014Muscle::updMuscleDynamicsInfo(
    const SimTK::State& s) const {

}

const MusclePotentialEnergyInfo&
Haeufle2014Muscle::getMusclePotentialEnergyInfo(
    const SimTK::State& s) const
{

}

MusclePotentialEnergyInfo& Haeufle2014Muscle::updMusclePotentialEnergyInfo(
    const SimTK::State& s) const {

}


double Haeufle2014Muscle::clampFiberLength(double lce) const {

}