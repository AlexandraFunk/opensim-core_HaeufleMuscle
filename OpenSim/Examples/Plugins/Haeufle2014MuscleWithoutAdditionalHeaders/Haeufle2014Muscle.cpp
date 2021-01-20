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
    constructProperty_exponent_descending_active_force_length(1.50);
    constructProperty_width_ascending_active_force_length(0.45);
    constructProperty_exponent_ascending_active_force_length(3.00);
    constructProperty_width_ascending_active_force_length(0.45);
    constructProperty_concentric_contraction_a_rel0(0.2);
    constructProperty_concentric_contraction_b_rel0(2.0);
    constructProperty_max_force_eccentric_extension(1.5);
    constructProperty_slopefactor(2.0);
    constructProperty_parallel_elastic_zero_length(0.95);
    constructProperty_parallel_elastic_exponent(2.5);
    constructProperty_parallel_elastic_force_rel_to_fmax(2.0);
    constructProperty_relative_stretch_at_nonlinear_linear_transition(0.0425);
    constructProperty_relative_stretch_at_linear_part(0.0170);
    constructProperty_force_at_nonlinear_linear_transition(getMaxIsometricForce() * 0.4); // dFsee0 = 0.4 * Fmax
    constructProperty_dse_damping_factor(0.3);
    constructProperty_rse_damping_factor(0.01);
    constructProperty_dpe_damping_factor(0);
    constructProperty_rpe_damping_factor(1);

    //TODO check if this is necessary?
    // setMinControl(get_minimum_activation());

}

void Haeufle2014Muscle::extendFinalizeFromProperties() 
{
    Super::extendFinalizeFromProperties();

    // Set the names of the muscle curves.
    const std::string& namePrefix = getName();

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
    OPENSIM_THROW_IF_FRMOBJ(get_concentric_contraction_a_rel0() <= 0,
            InvalidPropertyValue,
            getProperty_concentric_contraction_a_rel0().getName(),
            "The relative concentric contraction dynamics of CE parameter "
            "Arel0 must be greater than zero");
    OPENSIM_THROW_IF_FRMOBJ(get_concentric_contraction_b_rel0() <= 0,
            InvalidPropertyValue,
            getProperty_concentric_contraction_b_rel0().getName(),
            "The relative concentric contraction dynamics of CE parameter "
            "Brel0 must be greater than zero");
    OPENSIM_THROW_IF_FRMOBJ(get_max_force_eccentric_extension() <= 0,
            InvalidPropertyValue,
            getProperty_max_force_eccentric_extension().getName(),
            "The maximum eccentric extension force must be greater than zero");
    OPENSIM_THROW_IF_FRMOBJ(get_slopefactor() <= 0, InvalidPropertyValue,
            getProperty_slopefactor().getName(),
            "The corresponding slopefactor must be greater than zero");
    OPENSIM_THROW_IF_FRMOBJ(get_parallel_elastic_zero_length() <= 0,
            InvalidPropertyValue,
            getProperty_parallel_elastic_zero_length().getName(),
            "The zero length of the parallel elastic element must be greater "
            "than zero");
    OPENSIM_THROW_IF_FRMOBJ(get_parallel_elastic_exponent() <= 0,
            InvalidPropertyValue,
            getProperty_parallel_elastic_exponent().getName(),
            "The exponent of the parallel elastic element must be greater than "
            "zero");
    OPENSIM_THROW_IF_FRMOBJ(get_parallel_elastic_force_rel_to_fmax() <= 0,
            InvalidPropertyValue,
            getProperty_parallel_elastic_force_rel_to_fmax().getName(),
            "The elastic force relative to the maximum isometric force of the "
            "parallel elastic element must be greater than zero");
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
    OPENSIM_THROW_IF_FRMOBJ(get_dse_damping_factor() <= 0, InvalidPropertyValue,
            getProperty_dse_damping_factor().getName(),
            "The dse damping factor must be greater than zero")
    OPENSIM_THROW_IF_FRMOBJ(get_rse_damping_factor() <= 0, InvalidPropertyValue,
            getProperty_rse_damping_factor().getName(),
            "The rse damping factor must be greater than zero");

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

double Haeufle2014Muscle::getExponentDescendingActiveForceLength() const {
    return get_exponent_descending_active_force_length();
}

double Haeufle2014Muscle::getWidthDescendingActiveForceLength() const {
    return get_width_descending_active_force_length();
}

double Haeufle2014Muscle::getExponentAscendingActiveForceLength() const {
    return get_exponent_ascending_active_force_length();
}

double Haeufle2014Muscle::getWidthAscendingActiveForceLength() const {
    return get_width_ascending_active_force_length();
}

double Haeufle2014Muscle::getConcentricContractionARel0() const {
    return get_concentric_contraction_a_rel0();
}

double Haeufle2014Muscle::getConcentricContractionBRel0() const {
    return get_concentric_contraction_b_rel0();
}

double Haeufle2014Muscle::getMaxForceEccentricExtension() const {
    return get_max_force_eccentric_extension();
}

double Haeufle2014Muscle::getSlopeFactor() const {
    return get_slopefactor();
}

double Haeufle2014Muscle::getParallelElasticZeroLength() const {
    return get_parallel_elastic_zero_length();
}

double Haeufle2014Muscle::getParallelElasticExponent() const {
    return get_parallel_elastic_exponent();
}

double Haeufle2014Muscle::getParallelElasticForceRelToFmax() const {
    return get_parallel_elastic_force_rel_to_fmax();
}

double
Haeufle2014Muscle::getRelativeStretchAtNonlinearLinearTransition()
        const {
    return get_relative_stretch_at_nonlinear_linear_transition();
}

double Haeufle2014Muscle::getRelativeStretchAtLinearPart() const {
    return get_relative_stretch_at_linear_part();
}

double Haeufle2014Muscle::getForceAtNonlinearLinearTransition() const {
    return get_force_at_nonlinear_linear_transition();
}

double Haeufle2014Muscle::getDseDampingFactor() const {
    return get_dse_damping_factor();
}

double Haeufle2014Muscle::getRseDampingFactor() const {
    return get_rse_damping_factor();
}

    /** @returns The dse damping factor */
double Haeufle2014Muscle::getDpeDampingFactor() const {
    return get_dpe_damping_factor();
}

/** @returns The rse damping factor */
double Haeufle2014Muscle::getRpeDampingFactor() const {
    return get_rpe_damping_factor();
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

void Haeufle2014Muscle::setActiveForceLengthParameters(
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

void Haeufle2014Muscle::setFiberVelocityParameters(
        double aConcentricContractionARel0, double aConcentricContractionBRel0,
        double aMaxForceEccentricExtension, double aSlopeFactor) {
    set_concentric_contraction_a_rel0(aConcentricContractionARel0);
    set_concentric_contraction_b_rel0(aConcentricContractionBRel0);
    set_max_force_eccentric_extension(aMaxForceEccentricExtension);
    set_slopefactor(aSlopeFactor);
}

void Haeufle2014Muscle::setParallelElasticParameters(
        double aParallelElasticZeroLength, double aParallelElasticExponent,
        double aParallelElasticForceRelToFmax) {
    set_parallel_elastic_zero_length(aParallelElasticZeroLength);
    set_parallel_elastic_exponent(aParallelElasticExponent);
    set_parallel_elastic_force_rel_to_fmax(aParallelElasticForceRelToFmax);
}

void Haeufle2014Muscle::setForceAtNonlinearLinearTransition(
        double aForceAtNonlinearLinearTransition) {
    set_force_at_nonlinear_linear_transition(aForceAtNonlinearLinearTransition);
}

void Haeufle2014Muscle::
        setRelativeStretchAtNonlinearLinearTransition(
                double aRelativeStretchAtNonlinearLinearTransition) {
    set_relative_stretch_at_nonlinear_linear_transition(
            aRelativeStretchAtNonlinearLinearTransition);
}

void Haeufle2014Muscle::setRelativeStretchAtLinearPart(
        double aRelativeStretchAtLinearPart) 
{
    set_relative_stretch_at_linear_part(aRelativeStretchAtLinearPart);
}

void Haeufle2014Muscle::setTendonDampingParams(double aDseDampingFactor,
        double aRseDampingFactor) 
{
    set_dse_damping_factor(aDseDampingFactor);
    set_rse_damping_factor(aRseDampingFactor);
}

void Haeufle2014Muscle::setParallelDampingParams(
    double aDpeDampingFactor, double aRpeDampingFactor) 
{
    set_dpe_damping_factor(aDpeDampingFactor);
    set_rpe_damping_factor(aRpeDampingFactor);
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
        // Get the quantities that we've already computed.
        const MuscleLengthInfo& mli = getMuscleLengthInfo(s);

        double lceopt = getOptimalFiberLength();
        double gamma = getCalciumConcentration(s);
        double ldotMTC = getLengtheningSpeed(s);

        /** start calculating the activation for this model based on the
         * normalized ca concentration */
        double activation = getActivationModel().calculateActivation(
                gamma, mli.fiberLength);
        /** start calculating the values for the Haeufle Modell in the
         * following
         * order:
         * Fisom, Fpee, Fsee, Arel, Brel, D0, C2, C1, C0
         */
        double Fisom = calcFisom(mli.fiberLength);
        double Fpee = calcFpee(mli.fiberLength);
        double Fsee = calcFsee(mli.fiberLength);
        double Arel = calcArel(mli.fiberLength, activation, Fisom);
        double Brel = calcBrel(activation);
        /**
         * For calculating the coefficients C2, C1 and C0 use the pennated
         * Haeufle muscle description by Maria Hammer. Here the equation is
         * formulated the following:
         * 0 = ( C2 * (ldot/lceopt)^2 + C1 * (ldot/lceopt) + C0 ) / Fmax
         * first we solve for the normalized fiberLengthening (ldot/lceopt)
         * and finally divide this by the maximum isometric force of this
         * muscle Fmax
         */
        double C2woFmax = calcC2PenMaria(mli.fiberLength, mli.cosPennationAngle,
                activation, Fisom, Fpee);
        double C1woFmax = calcC1PenMaria(mli.fiberLength, ldotMTC,
                mli.cosPennationAngle, activation, Fisom, Fpee, Fsee);
        double C0woFmax = calcC0PenMaria(
                ldotMTC, mli.cosPennationAngle, activation, Fisom, Fpee, Fsee);

        // solve the equation 
        // TODO how to solve this?
        if (C2woFmax > 0) {
            double x1 = -C1woFmax +
                        sqrt(pow(C1woFmax, 2) - 4 * C2woFmax * C0woFmax) /
                                (2 * C2woFmax);
            double x2 = -C1woFmax -
                        sqrt(pow(C1woFmax, 2) - 4 * C2woFmax * C0woFmax) /
                                (2 * C2woFmax);
        }

        // excentric or concentric case? How to decide?
        double solution_placeholder = 0;

        // select the correct solution and divide it by Fmax to get the normalized Fiber velocity
        double normFiberVelocity =
                solution_placeholder / getMaxIsometricForce();
        double vce = normFiberVelocity * lceopt;

        double dphidt = getPennationModel().calcPennationAngularVelocity(
                tan(mli.pennationAngle), mli.fiberLength, vce);
        double dlceAT = getPennationModel().calcFiberVelocityAlongTendon(
                mli.fiberLength, vce, mli.sinPennationAngle,
                mli.cosPennationAngle, dphidt);
        double dmcldt = getLengtheningSpeed(s);
        double dtl = 0;

        if (!get_ignore_tendon_compliance()) {
            dtl = getPennationModel().calcTendonVelocity(mli.cosPennationAngle,
                    mli.sinPennationAngle, dphidt, mli.fiberLength, vce,
                    dmcldt);
        }

        // Populate the struct.
        fvi.fiberVelocity = vce;
        fvi.normFiberVelocity = normFiberVelocity;
        fvi.fiberVelocityAlongTendon = dlceAT;
        fvi.pennationAngularVelocity = dphidt;
        fvi.tendonVelocity = dtl;
        fvi.normTendonVelocity = dtl / getTendonSlackLength();
        // This one is not implemented in this model
        fvi.fiberForceVelocityMultiplier = SimTK::NaN;

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
        // Get the quantities that we've already computed.
        const MuscleLengthInfo& mli = getMuscleLengthInfo(s);
        const FiberVelocityInfo& fvi = getFiberVelocityInfo(s);

        double lceopt = getOptimalFiberLength();
        double gamma = getCalciumConcentration(s);
        double ldotMTC = getLengtheningSpeed(s);
        double Fmax = getMaxIsometricForce();

        /** start calculating the activation for this model based on the
         * normalized ca concentration */
        double activation = getActivationModel().calculateActivation(
                gamma, mli.fiberLength);

        double Fpee = calcFpee(mli.fiberLength);
        double Fsee = calcFsee(mli.fiberLength);
        double Fce = calcFceWithoutFmax(fvi.normFiberVelocity,
                             mli.normFiberLength, activation) *
                     Fmax;
        double FfiberAT = mli.cosPennationAngle * (Fce + Fpee);
        double Fsde = calcFsde(fvi.tendonVelocity, FfiberAT);
        double Fpde = calcFpde();
        


        mdi.activation = activation;
        mdi.fiberForce = Fce + Fpee;
        mdi.fiberForceAlongTendon = FfiberAT;
        mdi.normFiberForce = (Fce + Fpee) / Fmax;
        mdi.activeFiberForce = Fce;
        mdi.passiveFiberForce = Fpee;
        mdi.tendonForce = Fsee + Fsde;
        // mdi.normTendonForce = fse;

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

        // not implemente, is this necessary for this muscle model?

        // Implement TendonPotentialEnergy
        // integral über Kraft Längen Kurve bis aktuelle Tendon länge
        // nur in Pee und See

        // implementierung analog zu Millard

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
    // Fehlerausgabe weil Muskel schlackert wenn kleiner 0
    return max(lce, 0.0);
}


//==============================================================================
// HAEUFLE FORCE FUNCTIONS
//==============================================================================

double Haeufle2014Muscle::calcFisom(double FiberLength) const {
    // define variables and set them to values which will let the simulation
    // fail if they are not overwritten.
    double exponent_active_force_length = 0;
    double width_active_force_length = 0;
    // the optimal fiber length
    double lceopt = getOptimalFiberLength();
    double normFiberLength = FiberLength / lceopt;
    // if normalized fiber length is greater than 1 we are in the descending
    // domain
    if (normFiberLength > 1) {
        exponent_active_force_length = getExponentDescendingActiveForceLength();
        width_active_force_length = getWidthDescendingActiveForceLength();
    } else { // if normalized fiber length is smaller than 1 we are in the
             // ascending domain
        exponent_active_force_length = getExponentAscendingActiveForceLength();
        width_active_force_length = getWidthAscendingActiveForceLength();
    }
    return exp(-pow(abs((normFiberLength - 1) / width_active_force_length),
            exponent_active_force_length));
}

double Haeufle2014Muscle::calcFceWithoutFmax(double normFiberVelocity,
        double normFiberLength, double activation) const {
    double Fce = 0; // initial Fce to 0
    // check if this is a eccentric of concentric lengthening
    if (normFiberVelocity <= 0) { // concentric lengthening
        double Fisom = calcFisom(normFiberLength);
        double Arel = calcArel(normFiberLength, activation, Fisom);
        double Brel = calcBrel(activation);
        Fce = (activation * Fisom + Arel) / (1 - normFiberVelocity / Brel) -
              Arel;
    } else { // ecentric lengthening
        double Fisom = calcFisom(normFiberLength);
        double Arele = calcArele(activation, Fisom);
        double Brele = calcBrele(normFiberLength, activation, Fisom);
        Fce = (activation * Fisom + Arele) /
                      (1 - normFiberVelocity / Brele) -
              Arele;
    }
    return Fce;
}

double Haeufle2014Muscle::calcArel(double FiberLength,
        double activation, double Fisom) const {
    double Arel0 = getConcentricContractionARel0();
    double normFiberLength = FiberLength / getOptimalFiberLength();
    double Qarel = 1 / 4 * (1 + 3 * activation);
    double Larel = 0; // initialize to set total Arel to zero
    if (normFiberLength < 1) {
        Larel = 1;
    } else {
        Larel = Fisom;
    }
    return Arel0 * Qarel * Larel;
}

double Haeufle2014Muscle::calcBrel(double activation) const {
    double Brel0 = getConcentricContractionBRel0();
    double Qbrel = 1 / 7 * (3 + 4 * activation);
    // double Lbrel = 1; // not necessary
    return Brel0 * Qbrel;
}

double Haeufle2014Muscle::calcArele(
        double activation, double Fisom) const {
    double maxEccentricForce = getMaxForceEccentricExtension();
    double Arele = -maxEccentricForce * activation * Fisom;
    return Arele;
}

double Haeufle2014Muscle::calcBrele(double normFiberLength,
        double activation, double Fisom) const {
    double Arel = calcArel(
            normFiberLength, activation, Fisom);
    double Brel = calcBrel(activation);
    double slopefactor = getSlopeFactor();
    double maxEccentricForce = getMaxForceEccentricExtension();

    double Brele = Brel * (1 - maxEccentricForce) /
                   (slopefactor * (1 + Arel / (activation * Fisom)));
    return Brel;
}

double Haeufle2014Muscle::calcFpee(double FiberLength) const {
    double Lceopt = getOptimalFiberLength();
    double Lpee0 = getParallelElasticZeroLength() * Lceopt;
    double Fpee = 0; // standard case that FiberLength is smaller than Lpee0
    // check if Fiberlength is larger than or equal to Lpee0
    if (FiberLength >= Lpee0) {
        double Kpee = calcKPEE();
        Fpee = Kpee * pow((FiberLength - Lpee0), getParallelElasticExponent());
    }
    return Fpee;
}

double Haeufle2014Muscle::calcKPEE() const {
    double Lceopt = getOptimalFiberLength();
    double Fpee = getParallelElasticForceRelToFmax();
    double Lpee0 = getParallelElasticZeroLength();
    double nuepee = getParallelElasticExponent();
    double Fmax = getMaxIsometricForce();
    double Wdes = getWidthDescendingActiveForceLength();
    return (Fpee * Fmax / (pow(Lceopt * (Wdes + 1 - Lpee0), nuepee)));
}

double Haeufle2014Muscle::calcFsee(double aSerialElasticLength) const 
{
    double Lsee0 = getTendonSlackLength();
    double Fsee =
            0; // initialize Fsee for first case aSerialElasticLength < Lsee0
    double deltaUseenll = getRelativeStretchAtNonlinearLinearTransition();
    if (aSerialElasticLength > Lsee0 &&
            aSerialElasticLength <= (Lsee0 * (1 + deltaUseenll))) {
        double deltaFsee0 = getForceAtNonlinearLinearTransition();
        double nuesee = deltaUseenll / getRelativeStretchAtLinearPart();
        Fsee = deltaFsee0 *
               pow((aSerialElasticLength - Lsee0) / (deltaUseenll * Lsee0),
                       nuesee);
    } else if (aSerialElasticLength > (Lsee0 * (1 + deltaUseenll))) {
        double deltaFsee0 = getForceAtNonlinearLinearTransition();
        double deltaUseel = getRelativeStretchAtLinearPart();
        Fsee = deltaFsee0 *
               (1 + (aSerialElasticLength - Lsee0 * (1 + deltaUseenll)) /
                               (deltaUseel * Lsee0));
    }
    return Fsee;
}

double Haeufle2014Muscle::calcFsde(double serialElasticLengthVelocity,
        double totalFiberForceAlongTendon) const {

    double Dse = getDseDampingFactor();
    double Rse = getRseDampingFactor();
    double Fmax = getMaxIsometricForce();
    double lceopt = getOptimalFiberLength();
    double Arel0 = getConcentricContractionARel0();
    double Brel0 = getConcentricContractionBRel0();

    double Fsde = Dse * Fmax * Arel0 / (lceopt * Brel0) *
                  ((1 - Rse) * totalFiberForceAlongTendon / Fmax + Rse) *
                  serialElasticLengthVelocity;
    return Fsde;
}

double Haeufle2014Muscle::calcFpde(
    double lengthVelocity, double totalFiberForce) const 
{
    double Dpe = getDpeDampingFactor();
    double Rpe = getRpeDampingFactor();
    double Fmax = getMaxIsometricForce();
    double lceopt = getOptimalFiberLength();
    double Arel0 = getConcentricContractionARel0();
    double Brel0 = getConcentricContractionBRel0();

    double Fpde = Dpe * Fmax * Arel0 / (lceopt * Brel0) *
                  ((1 - Rpe) * totalFiberForce / Fmax + Rpe) * lengthVelocity;

    return Fpde;
}

double Haeufle2014Muscle::calcC2PenMaria(double fiberLength, double cosPenAngle,
        double activation, double Fisom, double Fpee) const {
    double Arel0 = getConcentricContractionARel0();
    double Brel0 = getConcentricContractionBRel0();
    double Arel = calcArel(fiberLength, activation, Fisom);
    double Brel = calcBrel(activation);
    double Dse = getDseDampingFactor();
    double Rse = getRseDampingFactor();
    double Fmax = getMaxIsometricForce();
    // necessary to implement these factors, since they are not demoa yet (Fpde)
    // ?
    double Dpe = 0;
    double Rpe = 0;
    double C2woFmax =
            Arel0 / (Brel0 * Brel) *
            (Dse * Rse / cosPenAngle + Dpe * Rpe -
                    (Dse - Dse * Rse + Dpe - Dpe * Rpe) * (Arel - Fpee / Fmax));
    return C2woFmax;
}

double Haeufle2014Muscle::calcC1PenMaria(double fiberLength, double ldotMTC,
        double cosPenAngle, double activation, double Fisom, double Fpee,
        double Fsee) const
{
    double Arel0 = getConcentricContractionARel0();
    double Brel0 = getConcentricContractionBRel0();
    double lceopt = getOptimalFiberLength();
    double Dse = getDseDampingFactor();
    double Rse = getRseDampingFactor();
    double Fmax = getMaxIsometricForce();
    double Arel = calcArel(fiberLength, activation, Fisom);
    double Brel = calcBrel(activation);
    // necessary to implement these factors, since they are not demoa yet (Fpde) ?
    double Dpe = 0;
    double Rpe = 0;
    double C1woFmax =
            Dse * Arel0 * ldotMTC / (Brel0 * Brel * lceopt) *
                    (cosPenAngle * (1 - Rse) * Arel - Rse) -
            Arel0 / Brel0 *
                    ((Dse - Dse * Rse + cosPenAngle * Dpe -
                             cosPenAngle * Dpe * Rpe) *
                                    activation * Fisom +
                            cosPenAngle * Dpe * Rpe + Rse * Dse / cosPenAngle) -
            cosPenAngle * Arel / Brel - Fsee / (Brel * Fmax) +
            cosPenAngle * Fpee / (Brel * Fmax);
    return C1woFmax;
}

double Haeufle2014Muscle::calcC0PenMaria(double ldotMTC, double cosPenAngle,
        double activation, double Fisom, double Fpee, double Fsee) const {
    double Arel0 = getConcentricContractionARel0();
    double Brel0 = getConcentricContractionBRel0();
    double lceopt = getOptimalFiberLength();
    double Dse = getDseDampingFactor();
    double Rse = getRseDampingFactor();
    double Fmax = getMaxIsometricForce();
    
    double C0woFmax =
            Dse * Arel0 * ldotMTC / (lceopt * Brel0) *
                    ((1 - Rse) * cosPenAngle * activation * Fisom + Rse) -
            cosPenAngle * activation * Fisom - cosPenAngle * Fpee / Fmax +
            Fsee / Fmax;
    return C0woFmax;
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
        double Fisom_lower = calcFisom(lower_border);
        double Fisom_middle = calcFisom(middle_border);
        // calculate Fpee with lower and middle_border 
        double Fpee_lower = calcFpee(lower_border);
        double Fpee_middle = calcFpee(middle_border);
        /** 
        * calculate Fsee with lower and middle_border
        * since lmtc = lsee + cos(alpha) * lce
        * -> lsee = lmtc - cos(alpha) * lce
        * For lower border this means: lce = lower_border
        * -> lsee = lmtc - cos(alpha) * lowerborder
        * FOr middle border this means: lce = middle_border
        * -> lsee = lmtc - cos(alpha) * lowerborder
        */
        double cosPenAng_lower =
                cos(getPennationModel().calcPennationAngle(lower_border));
        double lsee_lower = getPennationModel().calcTendonLength(
                cosPenAng_lower, lower_border, pathLength);
        double Fsee_lower = calcFsee(lsee_lower);
        double cosPenAng_middle =
                cos(getPennationModel().calcPennationAngle(middle_border));
        double lsee_middle = getPennationModel().calcTendonLength(
                cosPenAng_middle, middle_border, pathLength);
        double Fsee_middle = calcFsee(lsee_middle);

        // calculate the total lower and middle values:
        double initalEquilibrium_lower =
                cosPenAng_lower *
                        (aActivation * Fisom_lower * Fmax + Fpee_lower) -
                Fsee_lower;
        double initalEquilibrium_middle =
                cosPenAng_middle *
                        (aActivation * Fisom_middle * Fmax + Fpee_middle) -
                Fsee_middle;
        if ((initalEquilibrium_lower * initalEquilibrium_middle) < 0) {
            upper_border = middle_border;
        } else {
            lower_border = middle_border;
        }
        ferr = upper_border - lower_border;
    }

    // calculate final values:
    double lce_init = ((upper_border + lower_border) / 2);
    double cosPenAng_init =
            cos(getPennationModel().calcPennationAngle(lce_init));
    double lsee_init = getPennationModel().calcTendonLength(
            cosPenAng_init, lce_init, pathLength);
    double Fsee_init = calcFsee(lsee_init);

        // Populate the result map.
    ValuesFromInitMuscleState resultValues;

    if (abs(ferr) < aSolTolerance) { // The solution converged.

        resultValues["solution_error"] = ferr;
        resultValues["iterations"] = (double)iter;
        resultValues["fiber_length"] = lce_init;
        resultValues["tendon_force"] = Fsee_init;

        return std::pair<StatusFromInitMuscleState,
                ValuesFromInitMuscleState>(
                StatusFromInitMuscleState::Success_Converged,
                resultValues);
    }

    resultValues["solution_error"] = ferr;
    resultValues["iterations"] = (double)iter;
    resultValues["fiber_length"] = SimTK::NaN;
    resultValues["fiber_velocity"] = SimTK::NaN;
    resultValues["tendon_force"] = SimTK::NaN;

    return std::pair<StatusFromInitMuscleState, ValuesFromInitMuscleState>(
            StatusFromInitMuscleState::Failure_MaxIterationsReached,
            resultValues);
}
