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
    setNull();
    constructProperties();
    
    setName(aName);
    setMaxIsometricForce(aMaxIsometricForce);
    setOptimalFiberLength(aOptimalFiberLength);
    setTendonSlackLength(aTendonSlackLength);
    setPennationAngleAtOptimalFiberLength(aPennationAngle);
}

void Haeufle2014Muscle::setNull() 
{ 
    setAuthors("Maria Hammer, Mike Spahr"); 
}

void Haeufle2014Muscle::constructProperties() {
    constructProperty_default_calcium_concentration(0.05);
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

void Haeufle2014Muscle::extendFinalizeFromProperties() 
{
    Super::extendFinalizeFromProperties();

    // Set the names of the muscle curves.
    const std::string& namePrefix = getName();

    HaeufleActiveForceLengthCurve& falCurve = upd_HaeufleActiveForceLengthCurve();
    falCurve.setName(namePrefix + "_HaeufleActiveForceLengthCurve");

    HaeufleForceVelocityCurve& fvCurve = upd_HaeufleForceVelocityCurve();
    fvCurve.setName(namePrefix + "_HaeufleForceVelocityCurve");

    HaeufleFiberForceLengthCurve& fpeCurve = upd_HaeufleFiberForceLengthCurve();
    fpeCurve.setName(namePrefix + "_HaeufleFiberForceLengthCurve");

    HaeufleTendonForceLengthCurve& fseCurve = upd_HaeufleTendonForceLengthCurve();
    fseCurve.setName(namePrefix + "_HaeufleTendonForceLengthCurve");

    HaeufleTendonDampingCurve& fsdCurve = upd_HaeufleTendonDampingCurve();
    fsdCurve.setName(namePrefix + "_HaeufleTendonDampingCurve");

    // TODO add checks for too small fiber lenghts, etc... and throw errors when
    // claculations wouldns be possible!
    // throw errors if variables are in the wrong range 
    // check millard2012EquilibriumMuscle function for more info




    // Propagate properties down to pennation model subcomponent. If any of the
    // new property values are invalid, restore the subcomponent's current
    // property values (to avoid throwing again when the subcomponent's
    // extendFinalizeFromProperties() method is called directly) and then
    // re-throw the exception thrown by the subcomponent.
    auto& penMdl =
            updMemberSubcomponent<MuscleFixedWidthPennationModel>(penMdlIdx);
    MuscleFixedWidthPennationModel penMdlCopy(penMdl);
    penMdl.set_optimal_fiber_length(getOptimalFiberLength());
    penMdl.set_pennation_angle_at_optimal(
            getPennationAngleAtOptimalFiberLength());
    // this should be set automatically by the pennation model constructor
    // penMdl.set_maximum_pennation_angle(get_maximum_pennation_angle());
    try {
        penMdl.finalizeFromProperties();
    } catch (const InvalidPropertyValue&) {
        penMdl = penMdlCopy;
        throw;
    }

    // Propagate properties down to activation dynamics model subcomponent.
    // Handle invalid properties as above for pennation model.
    if (!get_ignore_activation_dynamics()) {
        auto& actMdl =
                updMemberSubcomponent<RockenfellerFirstOrderActivationDynamicModel>(
                        actMdlIdx);
        RockenfellerFirstOrderActivationDynamicModel actMdlCopy(actMdl);
        // set some properties of this model like this:
        // actMdl.set_
        try {
            actMdl.finalizeFromProperties();
        } catch (const InvalidPropertyValue&) {
            actMdl = actMdlCopy;
            throw;
        }
    }


    // TODO check if we need this for the first implementation of this model
    /**
    // Compute and store values that are used for clamping the fiber length.
    const double minActiveFiberLength =
            falCurve.getMinActiveFiberLength() * getOptimalFiberLength();
    const double minPennatedFiberLength = penMdl.getMinimumFiberLength();
    m_minimumFiberLength = max(SimTK::SignificantReal,
            max(minActiveFiberLength, minPennatedFiberLength));

    const double phi = penMdl.calcPennationAngle(m_minimumFiberLength);
    m_minimumFiberLengthAlongTendon =
            penMdl.calcFiberLengthAlongTendon(m_minimumFiberLength, cos(phi));
    */
}

//==============================================================================
// SCALING
//==============================================================================
void Haeufle2014Muscle::extendPostScale(
        const SimTK::State& s, const ScaleSet& scaleSet) {
    Super::extendPostScale(s, scaleSet);

    GeometryPath& path = upd_GeometryPath();
    if (path.getPreScaleLength(s) > 0.0) {
        double scaleFactor = path.getLength(s) / path.getPreScaleLength(s);
        upd_optimal_fiber_length() *= scaleFactor;
        upd_tendon_slack_length() *= scaleFactor;

        // Clear the pre-scale length that was stored in the GeometryPath.
        path.setPreScaleLength(s, 0.0);
    }
}


void Haeufle2014Muscle::extendConnectToModel(Model& model) 
{
    Super::extendConnectToModel(model);
}

void Haeufle2014Muscle::extendAddToSystem(
    SimTK::MultibodySystem& system) const 
{
    Super::extendAddToSystem(system);

    if (!get_ignore_activation_dynamics()) {
        addStateVariable(STATE_CALCIUM_CONCENTRATION);
    }
    if (!get_ignore_tendon_compliance()) {
        addStateVariable(STATE_FIBER_LENGTH_NAME);
    }
}

void Haeufle2014Muscle::extendInitStateFromProperties(SimTK::State& s) const 
{
    Super::extendInitStateFromProperties(s);

    if (!get_ignore_activation_dynamics()) {
        // sets the calcium concentration
        setActivation(s, getDefaultCalciumConcentration());
    }
    if (!get_ignore_tendon_compliance()) {
        setFiberLength(s, getDefaultFiberLength());
    }
}

void Haeufle2014Muscle::extendSetPropertiesFromState(const SimTK::State& s) 
{
    Super::extendSetPropertiesFromState(s);
    if (!get_ignore_activation_dynamics()) {
        setDefaultCalciumConcentration(
                getStateVariableValue(s, STATE_CALCIUM_CONCENTRATION));
    }
    if (!get_ignore_tendon_compliance()) {
        setDefaultFiberLength(
                getStateVariableValue(s, STATE_FIBER_LENGTH_NAME));
    }
}

void Haeufle2014Muscle::computeStateVariableDerivatives(
    const SimTK::State& s) const 
{
    // Activation dynamics if not ignored
    if (!get_ignore_activation_dynamics()) {
        double gammadot = 0;
        // if not disabled or overridden then compute its derivative
        if (appliesForce(s) && !isActuationOverridden(s)) {
            gammadot = getCalciumConcentrationDerivative(s);
        }
        setStateVariableDerivativeValue(s, STATE_CALCIUM_CONCENTRATION, gammadot);
    }

    // Fiber length is the next state (if it is a state at all)
    if (!get_ignore_tendon_compliance()) {
        double ldot = 0;
        // if not disabled or overridden then compute its derivative
        if (appliesForce(s) && !isActuationOverridden(s)) {
            ldot = getFiberVelocity(s);
        }
        setStateVariableDerivativeValue(s, STATE_FIBER_LENGTH_NAME, ldot);
    }

}

//==============================================================================
// GET METHODS
//==============================================================================

double Haeufle2014Muscle::getDefaultFiberLength() const 
{
    return get_default_fiber_length();
}

double Haeufle2014Muscle::getDefaultCalciumConcentration() const 
{
    return get_default_calcium_concentration();
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

double Haeufle2014Muscle::getCalciumConcentration(const SimTK::State& s) const {
    return getStateVariableValue(s, STATE_CALCIUM_CONCENTRATION);
}

double Haeufle2014Muscle::getCalciumConcentrationDerivative(const SimTK::State& s) const {
    if (get_ignore_activation_dynamics()) { return 0.0; }

    return getActivationModel().calcDerivative(
            getCalciumConcentration(s), getExcitation(s));
}

//==============================================================================
// SET METHODS
//==============================================================================

void Haeufle2014Muscle::setDefaultCalciumConcentration(
        double calciumConcentration) {
    set_default_calcium_concentration(calciumConcentration);
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
    SimTK::State& s, double activation) const 
{
    if (get_ignore_activation_dynamics()) {
        SimTK::Vector& controls(_model->updControls(s));
        setControls(SimTK::Vector(1, activation), controls);
        _model->setControls(s, controls);
    } else {
        setStateVariableValue(s, STATE_CALCIUM_CONCENTRATION,
                getActivationModel().clampGamma(activation));
    }
    markCacheVariableInvalid(s, _velInfoCV);
    markCacheVariableInvalid(s, _dynamicsInfoCV);
}

double Haeufle2014Muscle::computeActuation(const SimTK::State& s) const {
    const MuscleDynamicsInfo& mdi = getMuscleDynamicsInfo(s);
    setActuation(s, mdi.tendonForce);
    return mdi.tendonForce;
}

void Haeufle2014Muscle::computeInitialFiberEquilibrium(SimTK::State& s) const 
{
    // Initial activation and fiber length from input State, s.
    _model->getMultibodySystem().realize(s, SimTK::Stage::Velocity);
    double activation = getActivation(s);
    double pathLength = getLength(s);

    // Tolerance, in Newtons, of the desired equilibrium
    const double tol =
            max(1e-8 * getMaxIsometricForce(), SimTK::SignificantReal * 10);

    int maxIter = 200; // Should this be user settable?

    try {
        std::pair<StatusFromInitMuscleState, ValuesFromInitMuscleState> result =
                initMuscleState(activation, pathLength, tol, maxIter);

        switch (result.first) {

        case StatusFromInitMuscleState::Success_Converged:
            setActuation(s, result.second["tendon_force"]);
            setFiberLength(s, result.second["fiber_length"]);
            break;

        case StatusFromInitMuscleState::Warning_FiberAtLowerBound:
            log_warn("Haeufle2014Muscle static solution: '{}' is "
                     "at its minimum fiber length of {}.",
                    getName(), result.second["fiber_length"]);
            setActuation(s, result.second["tendon_force"]);
            setFiberLength(s, result.second["fiber_length"]);
            break;

        case StatusFromInitMuscleState::Failure_MaxIterationsReached:
            // Report internal variables and throw exception.
            std::ostringstream ss;
            ss << "\n  Solution error " << abs(result.second["solution_error"])
               << " exceeds tolerance of " << tol << "\n"
               << "  Newton iterations reached limit of " << maxIter << "\n"
               << "  Activation is " << activation << "\n"
               << "  Fiber length is " << result.second["fiber_length"] << "\n";
            OPENSIM_THROW_FRMOBJ(MuscleCannotEquilibrate, ss.str());
            break;
        }

    } catch (const std::exception& x) {
        OPENSIM_THROW_FRMOBJ(MuscleCannotEquilibrate,
                "Internal exception encountered.\n" + std::string{x.what()});
    }

}

void Haeufle2014Muscle::calcMuscleLengthInfo(
        const SimTK::State& s, MuscleLengthInfo& mli) const 
{
    try {
        // Get musculotendon actuator properties.
        double optFiberLength = getOptimalFiberLength();
        double tendonSlackLen = getTendonSlackLength();

        if (get_ignore_tendon_compliance()) { // rigid tendon
            mli.fiberLength =
                    clampFiberLength(getPennationModel().calcFiberLength(
                            getLength(s), tendonSlackLen));
        } else { // elastic tendon
            mli.fiberLength = clampFiberLength(
                    getStateVariableValue(s, STATE_FIBER_LENGTH_NAME));
        }

        mli.normFiberLength = mli.fiberLength / optFiberLength;
        mli.pennationAngle =
            getPennationModel().calcPennationAngle(mli.fiberLength);
        mli.cosPennationAngle = cos(mli.pennationAngle);
        mli.sinPennationAngle = sin(mli.pennationAngle);
        mli.fiberLengthAlongTendon = mli.fiberLength * mli.cosPennationAngle;

        // Necessary even for the rigid tendon, as it might have gone slack.
        mli.tendonLength = getPennationModel().calcTendonLength(
            mli.cosPennationAngle, mli.fiberLength, getLength(s));
        mli.normTendonLength = mli.tendonLength / tendonSlackLen;
        /** Tendon strain is defined using the elongation of the material
         * divided by its resting length.This is identical to the engineering
         * definition of strain.Thus a tendonStrain of 0.01 means that the
         * tendon is currently 1 % longer than its resting length.*/
        mli.tendonStrain = mli.normTendonLength - 1.0;

        // This model doesnt use this multipliers
        mli.fiberPassiveForceLengthMultiplier = SimTK::NaN;
        mli.fiberActiveForceLengthMultiplier = SimTK::NaN;

    }
    catch (const std::exception& x) {
        std::string msg = "Exception caught in Haeufle2014Muscle::"
            "calcMuscleLengthInfo from " +
            getName() + "\n" + x.what();
        throw OpenSim::Exception(msg);
    }
}

void Haeufle2014Muscle::calcFiberVelocityInfo(
    const SimTK::State& s, FiberVelocityInfo& fvi) const
{
    try {

    }
    catch (const std::exception& x) {
        std::string msg = "Exception caught in Haeufle2014Muscle::"
            "calcFiberVelocityInfo from " +
            getName() + "\n" + x.what();
        throw OpenSim::Exception(msg);
    }

}

void Haeufle2014Muscle::calcMuscleDynamicsInfo(
    const SimTK::State& s, MuscleDynamicsInfo& mdi) const
{
    try {

    }
    catch (const std::exception& x) {
        std::string msg = "Exception caught in Haeufle2014Muscle::"
            "calcMuscleDynamicsInfo from " +
            getName() + "\n" + x.what();
        throw OpenSim::Exception(msg);
    }
}

void Haeufle2014Muscle::calcMusclePotentialEnergyInfo(
    const SimTK::State& s, MusclePotentialEnergyInfo& mpei) const
{
    try {

    }
    catch (const std::exception& x) {
        std::string msg = "Exception caught in Haeufle2014Muscle::"
            "calcMusclePotentialEnergyInfo from " +
            getName() + "\n" + x.what();
        throw OpenSim::Exception(msg);
    }
}


double Haeufle2014Muscle::clampFiberLength(double lce) const
{
    // is this function necessary?
    return max(lce, 0.0);
}


std::pair<Haeufle2014Muscle::StatusFromInitMuscleState,
    Haeufle2014Muscle::ValuesFromInitMuscleState>
    Haeufle2014Muscle::initMuscleState(
        const double aActivation, const double pathLength, const double aSolTolerance,
        const int aMaxIterations) const
{
    /** implementing a simple Bisection method to solve the inital fiber length
     * problem:
     * cos(alpha) * ( q * Fisom * Fmax + Fpee) - Fsee = 0
     * within the interval (0, lMTC)
     */

     // get interval borders:
    double lower_border = 0;
    double upper_border = pathLength;

    // initialize iteration variable
    int iter = 0;
    double ferr = SimTK::MostPositiveReal; // Solution error

    while ((abs(ferr) > aSolTolerance) && (iter < aMaxIterations))
    {
        double middle_border = (lower_border + upper_border) / 2;

        double Fmax = getMaxIsometricForce();
        // calculate Fisom with lower and middle_border

        // calculate Fpee with lower and middle_border 

        // calculate Fsee with lower and middle_border
    }

}
