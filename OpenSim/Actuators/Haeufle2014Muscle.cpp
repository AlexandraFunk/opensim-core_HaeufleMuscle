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
    constructProperty_default_calcium_concentration(0.0);
    constructProperty_default_fiber_length(getOptimalFiberLength());
    constructProperty_fibre_exponent_descending_active_force_length(1.50);
    constructProperty_fibre_width_descending_active_force_length(0.45);
    constructProperty_fibre_exponent_ascending_active_force_length(3.00);
    constructProperty_fibre_width_ascending_active_force_length(0.45);
    constructProperty_fibre_active_force_velocity_arel0(0.2);
    constructProperty_fibre_active_force_velocity_brel0(2.0);
    constructProperty_fibre_max_eccentric_force_rel_to_fmax(1.5);
    constructProperty_fibre_eccentric_slopefactor(2.0);
    constructProperty_fibre_elastic_zero_length(0.95);
    constructProperty_fibre_elastic_exponent(2.5);
    constructProperty_fibre_elastic_force_rel_to_fmax(2.0);
    constructProperty_tendon_elastic_nonlinear_strain(0.0425);
    constructProperty_tendon_elastic_linear_strain(0.0170);
    constructProperty_tendon_elastic_force_rel_to_fmax(0.4); // before: dFsee0 = 0.4 * Fmax
    constructProperty_tendon_maximum_damping_dse(0.3);
    constructProperty_tendon_offset_damping_rse(0.01);
    constructProperty_fibre_maximum_damping_dpe(0.0);
    constructProperty_fibre_offset_damping_rpe(0.0);
    constructProperty_maximum_pennation_angle(acos(0.1)); // is this acos 0 or acos 0.1? 
    constructProperty_activation_Hatze_time_constant(11.3);
    constructProperty_activation_exponent(3.0);
    constructProperty_activation_optimal_calcium_concentration_fraction(
            5.27 * 1.37); // roh_0 * gamma_C
    constructProperty_activation_minimum(0.001); // Hatze constant
}

void Haeufle2014Muscle::extendFinalizeFromProperties() 
{
    Super::extendFinalizeFromProperties();

    // Set the names of the muscle curves.
    const std::string& namePrefix = getName();

    // Ensure property values are within appropriate ranges.
    OPENSIM_THROW_IF_FRMOBJ(get_fibre_exponent_descending_active_force_length() <= 0,
            InvalidPropertyValue,
            getProperty_fibre_exponent_descending_active_force_length().getName(),
            "The exponent of the descending branch of the active force length "
            "must be greater than zero");
    OPENSIM_THROW_IF_FRMOBJ(get_fibre_width_descending_active_force_length() <= 0,
            InvalidPropertyValue,
            getProperty_fibre_width_descending_active_force_length().getName(),
            "The width of the descending branch of the active force length "
            "curve must be greater than zero.");
    OPENSIM_THROW_IF_FRMOBJ(get_fibre_exponent_ascending_active_force_length() <= 0,
            InvalidPropertyValue,
            getProperty_fibre_exponent_ascending_active_force_length().getName(),
            "The exponent of the ascending branch of the active force length "
            "must be greater than zero.");
    OPENSIM_THROW_IF_FRMOBJ(get_fibre_width_ascending_active_force_length() <= 0,
            InvalidPropertyValue,
            getProperty_fibre_width_ascending_active_force_length().getName(),
            "The width of the ascending branch of the active force length must "
            "be greater than zero.");
    OPENSIM_THROW_IF_FRMOBJ(get_fibre_active_force_velocity_arel0() <= 0,
            InvalidPropertyValue,
            getProperty_fibre_active_force_velocity_arel0().getName(),
            "The relative concentric contraction dynamics of CE parameter "
            "Arel0 must be greater than zero");
    OPENSIM_THROW_IF_FRMOBJ(get_fibre_active_force_velocity_brel0() <= 0,
            InvalidPropertyValue,
            getProperty_fibre_active_force_velocity_brel0().getName(),
            "The relative concentric contraction dynamics of CE parameter "
            "Brel0 must be greater than zero");
    OPENSIM_THROW_IF_FRMOBJ(get_fibre_max_eccentric_force_rel_to_fmax() <= 0,
            InvalidPropertyValue,
            getProperty_fibre_max_eccentric_force_rel_to_fmax().getName(),
            "The maximum eccentric extension force must be greater than zero");
    OPENSIM_THROW_IF_FRMOBJ(get_fibre_eccentric_slopefactor() <= 0, InvalidPropertyValue,
            getProperty_fibre_eccentric_slopefactor().getName(),
            "The corresponding slopefactor must be greater than zero");
    OPENSIM_THROW_IF_FRMOBJ(get_fibre_elastic_zero_length() <= 0,
            InvalidPropertyValue,
            getProperty_fibre_elastic_zero_length().getName(),
            "The zero length of the parallel elastic element must be greater "
            "than zero");
    OPENSIM_THROW_IF_FRMOBJ(get_fibre_elastic_exponent() <= 0,
            InvalidPropertyValue,
            getProperty_fibre_elastic_exponent().getName(),
            "The exponent of the parallel elastic element must be greater than "
            "zero");
    OPENSIM_THROW_IF_FRMOBJ(get_fibre_elastic_force_rel_to_fmax() <= 0,
            InvalidPropertyValue,
            getProperty_fibre_elastic_force_rel_to_fmax().getName(),
            "The elastic force relative to the maximum isometric force of the "
            "parallel elastic element must be greater than zero");
    OPENSIM_THROW_IF_FRMOBJ(get_tendon_slack_length() <= 0,
            InvalidPropertyValue, getProperty_tendon_slack_length().getName(),
            "The serial elastic rest length must be greater than zero");
    OPENSIM_THROW_IF_FRMOBJ(
            get_tendon_elastic_nonlinear_strain() <= 0,
            InvalidPropertyValue,
            getProperty_tendon_elastic_nonlinear_strain()
                    .getName(),
            "The relative stretch at nonlinear/linear transition must be "
            "greater than zero");
    OPENSIM_THROW_IF_FRMOBJ(get_tendon_elastic_linear_strain() <= 0,
            InvalidPropertyValue,
            getProperty_tendon_elastic_linear_strain().getName(),
            "The relative stretch in the linear part must be greater than "
            "zero");
    OPENSIM_THROW_IF_FRMOBJ(get_tendon_elastic_force_rel_to_fmax() <= 0,
            InvalidPropertyValue,
            getProperty_tendon_elastic_force_rel_to_fmax().getName(),
            "The force at the nonlinear/linear transition must be greater than "
            "zero");
    OPENSIM_THROW_IF_FRMOBJ(get_tendon_maximum_damping_dse() < 0, InvalidPropertyValue,
            getProperty_tendon_maximum_damping_dse().getName(),
            "The dse damping factor must be greater than zero")
    OPENSIM_THROW_IF_FRMOBJ(get_tendon_offset_damping_rse() < 0, InvalidPropertyValue,
            getProperty_tendon_offset_damping_rse().getName(),
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
    penMdl.set_maximum_pennation_angle(get_maximum_pennation_angle());
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
        actMdl.set_minimum_gamma(get_default_calcium_concentration());
        actMdl.set_optimal_fiber_length(get_optimal_fiber_length());
        actMdl.set_activation_minimum(get_activation_minimum());
        actMdl.set_activation_exponent(get_activation_exponent());
        actMdl.set_activation_optimal_calcium_concentration_fraction(
                get_activation_optimal_calcium_concentration_fraction());
        actMdl.set_activation_Hatze_time_constant(
                get_activation_Hatze_time_constant());
        try {
            actMdl.finalizeFromProperties();
        } catch (const InvalidPropertyValue&) {
            actMdl = actMdlCopy;
            throw;
        }
    }

    // set maximum contraction velocity
    double vmax = calcVmax();
    setMaxContractionVelocity(vmax);

    // what shall be the minimum fiber length from penMdl??
    m_minimumFiberLength = 0;

    const double phi = penMdl.calcPennationAngle(m_minimumFiberLength);
    m_minimumFiberLengthAlongTendon =
            penMdl.calcFiberLengthAlongTendon(m_minimumFiberLength, cos(phi));
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
        double gammadot = 0.0;
        // if not disabled or overridden then compute its derivative
        if (appliesForce(s) && !isActuationOverridden(s)) {
            gammadot = getCalciumConcentrationDerivative(s);
        }
        setStateVariableDerivativeValue(s, STATE_CALCIUM_CONCENTRATION, gammadot);
    }

    // Fiber length is the next state (if it is a state at all)
    if (!get_ignore_tendon_compliance()) {
        double ldot = 0.0;
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

double Haeufle2014Muscle::getFibreExponentDescendingActiveForceLength() const {
    return get_fibre_exponent_descending_active_force_length();
}

double Haeufle2014Muscle::getFibreWidthDescendingActiveForceLength() const {
    return get_fibre_width_descending_active_force_length();
}

double Haeufle2014Muscle::getFibreExponentAscendingActiveForceLength() const {
    return get_fibre_exponent_ascending_active_force_length();
}

double Haeufle2014Muscle::getFibreWidthAscendingActiveForceLength() const {
    return get_fibre_width_ascending_active_force_length();
}

double Haeufle2014Muscle::getFibreActiveForceVelocityArel0() const {
    return get_fibre_active_force_velocity_arel0();
}

double Haeufle2014Muscle::getFibreActiveForceVelocityBrel0() const {
    return get_fibre_active_force_velocity_brel0();
}

double Haeufle2014Muscle::getFibreMaxEccentricForceRelToFmax() const {
    return get_fibre_max_eccentric_force_rel_to_fmax();
}

double Haeufle2014Muscle::getFibreEccentricSlopefactor() const {
    return get_fibre_eccentric_slopefactor();
}

double Haeufle2014Muscle::getFibreElasticZeroLength() const {
    return get_fibre_elastic_zero_length();
}

double Haeufle2014Muscle::getFibreElasticExponent() const {
    return get_fibre_elastic_exponent();
}

double Haeufle2014Muscle::getFibreElasticForceRelToFmax() const {
    return get_fibre_elastic_force_rel_to_fmax();
}

double
Haeufle2014Muscle::getTendonElasticNonlinearStrain()
        const {
    return get_tendon_elastic_nonlinear_strain();
}

double Haeufle2014Muscle::getTendonElasticLinearStrain() const {
    return get_tendon_elastic_linear_strain();
}

double Haeufle2014Muscle::getTendonElasticForceRelToFmax() const {
    return get_tendon_elastic_force_rel_to_fmax();
}

double Haeufle2014Muscle::getTendonMaximumDampingDse() const {
    return get_tendon_maximum_damping_dse();
}

double Haeufle2014Muscle::getTendonOffsetDampingRse() const {
    return get_tendon_offset_damping_rse();
}

double Haeufle2014Muscle::getFibreMaximumDampingDpe() const {
    return get_fibre_maximum_damping_dpe();
}

double Haeufle2014Muscle::getFibreOffsetDampingRpe() const {
    return get_fibre_offset_damping_rpe();
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

double Haeufle2014Muscle::getMinimumFiberLength() const {
    return m_minimumFiberLength;
}

double Haeufle2014Muscle::getMinimumFiberLengthAlongTendon() const {
    return m_minimumFiberLengthAlongTendon;
}

double Haeufle2014Muscle::getCalciumConcentrationDerivative(const SimTK::State& s) const 
{
    if (get_ignore_activation_dynamics()) { return 0.0; }

    return getActivationModel().calcDerivative(
            getCalciumConcentration(s), getExcitation(s));
}

double Haeufle2014Muscle::getFpee(const SimTK::State& s) const {
    return getMuscleDynamicsInfo(s).userDefinedDynamicsExtras[0];
}

double Haeufle2014Muscle::getFpde(const SimTK::State& s) const {
    return getMuscleDynamicsInfo(s).userDefinedDynamicsExtras[1];
}

double Haeufle2014Muscle::getFsee(const SimTK::State& s) const {
    return getMuscleDynamicsInfo(s).userDefinedDynamicsExtras[2];
}

double Haeufle2014Muscle::getFsde(const SimTK::State& s) const {
    return getMuscleDynamicsInfo(s).userDefinedDynamicsExtras[3];
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
        double fibreExponentDescendingActiveForceLength,
        double fibreWidthDescendingActiveForceLength,
        double fibreExponentAscendingActiveForceLength,
        double fibreWidthAscendingActiveForceLength) {
    set_fibre_exponent_descending_active_force_length(
            fibreExponentDescendingActiveForceLength);
    set_fibre_width_descending_active_force_length(fibreWidthDescendingActiveForceLength);
    set_fibre_exponent_ascending_active_force_length(
            fibreExponentAscendingActiveForceLength);
    set_fibre_width_ascending_active_force_length(fibreWidthAscendingActiveForceLength);
}

void Haeufle2014Muscle::setFiberVelocityParameters(
        double aFibreActiveForceVelocityArel0, double aFibreActiveForceVelocityBrel0,
        double aFibreMaxEccentricForceRelToFmax,
        double aFibreEccentricSlopefactor) {
    set_fibre_active_force_velocity_arel0(aFibreActiveForceVelocityArel0);
    set_fibre_active_force_velocity_brel0(aFibreActiveForceVelocityBrel0);
    set_fibre_max_eccentric_force_rel_to_fmax(aFibreMaxEccentricForceRelToFmax);
    set_fibre_eccentric_slopefactor(aFibreEccentricSlopefactor);
}

void Haeufle2014Muscle::setParallelElasticParameters(
        double aFibreElasticZeroLength, double aFibreElasticExponent,
        double aFibreElasticForceRelToFmax) {
    set_fibre_elastic_zero_length(aFibreElasticZeroLength);
    set_fibre_elastic_exponent(aFibreElasticExponent);
    set_fibre_elastic_force_rel_to_fmax(aFibreElasticForceRelToFmax);
}

void Haeufle2014Muscle::setTendonElasticForceRelToFmax(
        double aTendonElasticForceRelToFmax) {
    set_tendon_elastic_force_rel_to_fmax(
            aTendonElasticForceRelToFmax);
}

void Haeufle2014Muscle::
        setTendonElasticNonlinearStrain(
                double aTendonElasticNonlinearStrain) {
    set_tendon_elastic_nonlinear_strain(
            aTendonElasticNonlinearStrain);
}

void Haeufle2014Muscle::setTendonElasticLinearStrain(
        double aTendonElasticLinearStrain) 
{
    set_tendon_elastic_linear_strain(aTendonElasticLinearStrain);
}

void Haeufle2014Muscle::setTendonDampingParams(double aTendonMaximumDampingDse,
        double aTendonOffsetDampingRse) 
{
    set_tendon_maximum_damping_dse(aTendonMaximumDampingDse);
    set_tendon_offset_damping_rse(aTendonOffsetDampingRse);
}

void Haeufle2014Muscle::setParallelDampingParams(
    double aFibreMaximumDampingDpe, double aFibreOffsetDampingRpe) 
{
    set_fibre_maximum_damping_dpe(aFibreMaximumDampingDpe);
    set_fibre_offset_damping_rpe(aFibreOffsetDampingRpe);
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
// is necessary since it is implemented in the Muscle.h interface but in this case
// it doesnt set the activation rather then the calcium concentration stat variable
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

    double pathLength = getLength(s);

    // get initial excitation
    double excitation = getExcitation(s);

    // Tolerance, in Newtons, of the desired equilibrium
    const double tol = 1e-9;

    int maxIter = 200; // Should this be user settable?

    try {
        std::pair<StatusFromInitMuscleState, ValuesFromInitMuscleState> result =
                initMuscleState(pathLength, excitation, tol, maxIter);
        switch (result.first) {

        case StatusFromInitMuscleState::Success_Converged:
            setActuation(s, result.second["tendon_force"]);
            setFiberLength(s, result.second["fiber_length"]);
            setActivation(s, excitation);
            break;
        case StatusFromInitMuscleState::Failure_MaxIterationsReached:
            // Report internal variables and throw exception.
            std::ostringstream ss;
            ss << "\n  Solution error " << abs(result.second["solution_error"])
               << " exceeds tolerance of " << tol << "\n"
               << "  Newton iterations reached limit of " << maxIter << "\n"
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
        mli.fiberLengthAlongTendon =
                getPennationModel().calcFiberLengthAlongTendon(
                        mli.fiberLength, mli.cosPennationAngle);

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
        double activation = getActivationModel().clampActivation(
                getActivationModel().calculateActivation(
                        gamma, mli.fiberLength));
        /** start calculating the values for the Haeufle Modell in the
         * following
         * order:
         * Fisom, Fpee, Fsee, Arel, Brel, D0, C2, C1, C0
         */
        double Fisom = calcFisom(mli.fiberLength);
        double Fpee = calcFpee(mli.fiberLength);
        double Fsee = calcFsee(mli.tendonLength);

        /**
         * For calculating the coefficients C2, C1 and C0 use the pennated
         * Haeufle muscle description by Maria Hammer. Here the equation is
         * formulated the following:
         * 0 = ( C2 * (ldot/lceopt)^2 + C1 * (ldot/lceopt) + C0 ) / Fmax
         * first we solve for the normalized fiberLengthening (ldot/lceopt)
         * and finally divide this by the maximum isometric force of this
         * muscle Fmax
         */

        bool uniqueSolution = false;
        bool concentric_case;
        int error_case = 0;
        double lcedot = 0.0;
        double Arel = getFibreActiveForceVelocityArel0();
        double Brel = getFibreActiveForceVelocityBrel0();

        // ldotMTC <= 0 probably concentric movement (start calculating
        // with Arel)
        if (ldotMTC <= 0) {
            concentric_case = true;
        } else {
            concentric_case = false;
        }

        for (int loopvar = 0; loopvar < 2; loopvar++) {
            if (concentric_case) {
                Arel = calcArel(mli.fiberLength, activation, Fisom);
                Brel = calcBrel(activation);
            } else {
                // eccentric case
                Arel = calcArele(activation, Fisom);
                Brel = calcBrele(mli.fiberLength, activation, Fisom);
            }

            // start calculating C2, C1, C0 and D
            double C2dashed = calcC2dash(
                    mli.cosPennationAngle, activation, Fpee, Brel, Arel);
            double C1dashed = calcC1dash(ldotMTC, mli.cosPennationAngle,
                    activation, Fisom, Fpee, Fsee, Arel, Brel);
            double C0dashed = calcC0dash(ldotMTC, mli.cosPennationAngle,
                    activation, Fisom, Fpee, Fsee);
            double Ddashed = pow(C1dashed, 2.0) - 4.0 * C0dashed * C2dashed;

            if (abs(C2dashed) < SimTK::SignificantReal) {
                if (abs(C1dashed) < SimTK::SignificantReal) {
                    if (loopvar == 0) {
                        error_case = 1;
                        concentric_case = (!concentric_case);
                        continue;
                    } else {
                        break;
                    }
                } else {
                    lcedot = -lceopt * C0dashed / C1dashed;
                }
            } else {
                if (Ddashed < -SimTK::SignificantReal) { // this would lead to
                                                         // complex lcedot
                    if (loopvar == 0) {
                        error_case = 2;
                        concentric_case = (!concentric_case);
                        continue;
                    } else {
                        break;
                    }
                } else { // solve quadratic equation for lcedot
                    lcedot = lceopt * (-C1dashed - sqrt(Ddashed)) /
                             (2.0 * C2dashed);
                }
            }
            if (((concentric_case == true) &&
                        (lcedot <= SimTK::SignificantReal)) ||
                    ((concentric_case == false) &&
                            (lcedot >= -SimTK::SignificantReal))) {
                error_case = 0;
                break;
            } else {
                if (loopvar == 0) {
                    error_case = 3;
                    concentric_case = (!concentric_case);
                }
            }
        }

        if (error_case == 0) {
            // nothing to do here since everything is working fine
        } else if (error_case == 1) {
            log_warn("Haeufle2014Muscle::calcFiberVelocityInfo(): no solution "
                     "found for lcedot. C2=C1=0. Set lcedot=0.");
            lcedot = 0.0;
        } else if (error_case == 2) {
            log_warn("Haeufle2014Muscle::calcFiberVelocityInfo(): complex "
                     "solution found for lcedot. Set lcedot=0.");
            lcedot = 0;
        } else if (error_case == 3) {
            log_warn(
                    "Haeufle2014Muscle::calcFiberVelocityInfo(): no concentric "
                    "or eccentric solution found for lcedot. Set lcedot=0.");
            lcedot = 0.0;
        } else {
            log_warn("Haeufle2014Muscle::calcFiberVelocityInfo(): something is "
                     "wrong, error_case is out of range");
            lcedot = 0.0;
        }

        // get the seldom case that if Fce < 0 or Fce < Fisom(lce)*a0*Fmax
        // calculation of lcedot changes
        double Fmax = getMaxIsometricForce();
        double Fce = calcNormFce(lcedot, mli.fiberLength, activation) * Fmax;
        double minActivation = getActivationModel().get_activation_minimum();
        if (Fce < 0 || Fce < Fisom * minActivation * Fmax) {
            double Ffibre = Fpee / Fmax + minActivation * Fisom;
            double vmax = getMaxContractionVelocity();
            double Dse = getTendonMaximumDampingDse();
            double Rse = getTendonOffsetDampingRse();
            double Dpe = getFibreMaximumDampingDpe();
            double Rpe = getFibreOffsetDampingRpe();
            double cosAlpha = mli.cosPennationAngle;
            lcedot = (-cosAlpha * Ffibre * vmax + vmax * Fsee / Fmax +
                             ldotMTC * Dse *
                                     ((1 - Rse) * cosAlpha * Ffibre + Rse)) /
                     (Dse * ((1 - Rse) * Ffibre + Rse / cosAlpha) +
                             Dpe * cosAlpha * ((1 - Rpe) * Fpee / Fmax + Rpe));
        }


        double normFiberVelocity = lcedot / getMaxContractionVelocity();

        double dphidt = getPennationModel().calcPennationAngularVelocity(
                tan(mli.pennationAngle), mli.fiberLength, lcedot);
        double dlceAT = getPennationModel().calcFiberVelocityAlongTendon(
                mli.fiberLength, lcedot, mli.sinPennationAngle,
                mli.cosPennationAngle, dphidt);
        double dmcldt = getLengtheningSpeed(s);
        double dtl = 0;

        if (!get_ignore_tendon_compliance()) {
            dtl = getPennationModel().calcTendonVelocity(mli.cosPennationAngle,
                    mli.sinPennationAngle, dphidt, mli.fiberLength, lcedot,
                    dmcldt);
        }

        // Populate the struct.
        fvi.fiberVelocity = lcedot;
        fvi.fiberVelocityAlongTendon = dlceAT;
        fvi.pennationAngularVelocity = dphidt;
        fvi.tendonVelocity = dtl;
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
        double activation = getActivationModel().clampActivation(
                getActivationModel().calculateActivation(
                        gamma, mli.fiberLength));

        double Fpee = calcFpee(mli.fiberLength);
        double Fsee = calcFsee(mli.tendonLength);
        double Fce = calcNormFce(fvi.fiberVelocity,
                             mli.fiberLength, activation) *
                     Fmax;
        // get the seldom case that if Fce < 0 or Fce < Fisom(lce)*a0*Fmax
        // calculation of lcedot changes
        double Fisom = calcFisom(mli.fiberLength);
        double minActivation = getActivationModel().get_activation_minimum();
        if (Fce < 0 || Fce < Fisom * minActivation * Fmax) {
            Fce = Fisom * minActivation * Fmax;
        }

        double Fsde = calcFsde(fvi.tendonVelocity, mli.cosPennationAngle, Fpee, Fce);
        double Fpde = calcFpde(fvi.fiberVelocity, Fpee);

        mdi.activation = activation;
        mdi.fiberForce = Fce + Fpee + Fpde;
        mdi.fiberForceAlongTendon = mli.cosPennationAngle * (Fce + Fpee + Fpde);
        mdi.normFiberForce = (Fce + Fpee + Fpde) / Fmax;
        mdi.activeFiberForce = Fce;
        mdi.passiveFiberForce = Fpee + Fpde;
        mdi.tendonForce = Fsee + Fsde;
        mdi.normTendonForce = mdi.tendonForce / Fmax;

        // Set values which are not implemented in this model to NaN
        mdi.fiberStiffness = SimTK::NaN;
        mdi.fiberStiffnessAlongTendon = SimTK::NaN; 
        mdi.tendonStiffness = SimTK::NaN;
        mdi.muscleStiffness = SimTK::NaN; 
        mdi.fiberActivePower = SimTK::NaN;
        mdi.fiberPassivePower = SimTK::NaN; 
        mdi.tendonPower = SimTK::NaN;              
        mdi.musclePower = SimTK::NaN;               

        // to make all forces accessible
        mdi.userDefinedDynamicsExtras.resize(4);
        mdi.userDefinedDynamicsExtras[0] = Fpee;
        mdi.userDefinedDynamicsExtras[1] = Fpde;
        mdi.userDefinedDynamicsExtras[2] = Fsee;
        mdi.userDefinedDynamicsExtras[3] = Fsde;
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
        // Get the quantities that we've already computed.
        const MuscleLengthInfo& mli = getMuscleLengthInfo(s);

        // calculate integral for Fpee depending on current fiber length
        mpei.fiberPotentialEnergy = calcIntegralFpee(mli.fiberLength);
        // calculate integral for Fsee depending on current tendon length
        mpei.tendonPotentialEnergy = calcIntegralFsee(mli.tendonLength);
        // muscle potential enery is equal to the sum of the 
        // potential energy of the fiber and the tendon
        mpei.musclePotentialEnergy =
                mpei.fiberPotentialEnergy + mpei.tendonPotentialEnergy;
    }
    catch (const std::exception& x) {
        std::string msg = "Exception caught in Haeufle2014Muscle::"
            "calcMusclePotentialEnergyInfo from " +
            getName() + "\n" + x.what();
        throw OpenSim::Exception(msg);
    }
}

double Haeufle2014Muscle::calcVmax() const 
{
    double Arel0 = getFibreActiveForceVelocityArel0();
    double Brel0 = getFibreActiveForceVelocityBrel0();
    double lceopt = getOptimalFiberLength();

    double vmax = Brel0 * lceopt / Arel0;

    return vmax;
}

double Haeufle2014Muscle::clampFiberLength(double lce) const
{
    if (lce < getMinimumFiberLength()) {
        log_warn("Exception caught in Haeufle2014Muscle:: clampFiberLength "
                 "from {} \n Fiber length can't be smaller than 0",
                getName());
    }
    return max(lce, getMinimumFiberLength());
}


//==============================================================================
// HAEUFLE FORCE FUNCTIONS
//==============================================================================

double Haeufle2014Muscle::calcFisom(double fiberLength) const {
    // define variables and set them to values which will let the simulation
    // fail if they are not overwritten.
    double exponent_active_force_length = 0.0;
    double width_active_force_length = 0.0;
    // the optimal fiber length
    double lceopt = getOptimalFiberLength();
    double normFiberLength = fiberLength / lceopt;
    // if normalized fiber length is greater than 1 we are in the descending
    // domain
    if (normFiberLength > 1) {
        exponent_active_force_length = getFibreExponentDescendingActiveForceLength();
        width_active_force_length = getFibreWidthDescendingActiveForceLength();
    } else { // if normalized fiber length is smaller than 1 we are in the
             // ascending domain
        exponent_active_force_length = getFibreExponentAscendingActiveForceLength();
        width_active_force_length = getFibreWidthAscendingActiveForceLength();
    }
    return exp(-pow(abs((normFiberLength - 1.0) / width_active_force_length),
            exponent_active_force_length));
}

double Haeufle2014Muscle::calcNormFce(double fiberVelocity,
        double fiberLength, double activation) const {
    double Fce = 0; // initial Fce to 0
    double lceopt = getOptimalFiberLength();
    // check if this is a eccentric of concentric lengthening
    if (fiberVelocity <= 0) { // concentric lengthening
        double Fisom = calcFisom(fiberLength);
        double Arel = calcArel(fiberLength, activation, Fisom);
        double Brel = calcBrel(activation);
        Fce = (activation * Fisom + Arel) /
                      (1.0 - fiberVelocity / (Brel * lceopt)) -
              Arel;
    } else { // ecentric lengthening
        double Fisom = calcFisom(fiberLength);
        double Arele = calcArele(activation, Fisom);
        double Brele = calcBrele(fiberLength, activation, Fisom);
        Fce = (activation * Fisom + Arele) /
                      (1.0 - fiberVelocity / (Brele * lceopt)) -
              Arele;
    }
    return Fce;
}

double Haeufle2014Muscle::calcArel(double fiberLength,
        double activation, double Fisom) const {
    double Arel0 = getFibreActiveForceVelocityArel0();
    double normFiberLength = fiberLength / getOptimalFiberLength();
    double Qarel = 1.0 / 4.0 * (1.0 + 3.0 * activation);
    double Larel = (normFiberLength < 1) ? 1.0 : Fisom;
    double Arel = Arel0 * Qarel * Larel;
    return Arel;
}

double Haeufle2014Muscle::calcBrel(double activation) const {
    double Brel0 = getFibreActiveForceVelocityBrel0();
    double Qbrel = 1.0 / 7.0 * (3.0 + 4.0 * activation);
    // double Lbrel = 1; // not necessary due to multiplication
    double Brel = Brel0 * Qbrel;
    return Brel;
}

double Haeufle2014Muscle::calcArele(
        double activation, double Fisom) const {
    double maxEccentricForce = getFibreMaxEccentricForceRelToFmax();
    double Arele = -maxEccentricForce * activation * Fisom;
    return Arele;
}

double Haeufle2014Muscle::calcBrele(double fiberLength,
        double activation, double Fisom) const {
    double Arel = calcArel(
            fiberLength, activation, Fisom);
    double Brel = calcBrel(activation);
    double slopefactor = getFibreEccentricSlopefactor();
    double maxEccentricForce = getFibreMaxEccentricForceRelToFmax();
    double Brele = Brel * (1.0 - maxEccentricForce) /
                   (slopefactor * (1.0 + Arel / (activation * Fisom)));
    return Brele;
}

double Haeufle2014Muscle::calcFpee(double fiberLength) const {
    double Lceopt = getOptimalFiberLength();
    double Lpee0 = getFibreElasticZeroLength();
    double Fpee = 0.0; // standard case that FiberLength is smaller than Lpee0 * lceopt
    if (fiberLength >= (Lpee0 * Lceopt)) {
        // double Kpee = calcKPEE(); // deprecated
        double deltaFpee =
                getFibreElasticForceRelToFmax() * getMaxIsometricForce();
        double nuepee = getFibreElasticExponent();
        double Wdes = getFibreWidthDescendingActiveForceLength();
        Fpee = deltaFpee *
               pow((fiberLength / Lceopt - Lpee0) / (Wdes + 1 - Lpee0), nuepee);
    }
    return Fpee;
}

double Haeufle2014Muscle::calcIntegralFpee(double fiberLength) const 
{
    double Lceopt = getOptimalFiberLength();
    double Lpee0 = getFibreElasticZeroLength();
    double IntegralFpee = 0.0; // standard case that FiberLength is smaller than Lpee0
    if (fiberLength >= (Lpee0*Lceopt)) 
    {
        double Fpee = calcFpee(fiberLength);
        double nuepee = getFibreElasticExponent();
        IntegralFpee = Fpee * (fiberLength - Lpee0*Lceopt) / (nuepee + 1.0);      
    }
    return IntegralFpee;
}

double Haeufle2014Muscle::calcFsee(double aSerialElasticLength) const 
{
    double Lsee0 = getTendonSlackLength();
    // initialize Fsee for first case aSerialElasticLength < Lsee0
    double Fsee = 0.0; 
    double deltaUseenll = getTendonElasticNonlinearStrain();
    if (aSerialElasticLength > Lsee0 &&
            aSerialElasticLength <= (Lsee0 * (1.0 + deltaUseenll))) {
        double deltaFsee0 = getTendonElasticForceRelToFmax() * getMaxIsometricForce();
        double nuesee = deltaUseenll / getTendonElasticLinearStrain();
        Fsee = deltaFsee0 *
               pow((aSerialElasticLength - Lsee0) / (deltaUseenll * Lsee0),
                       nuesee);
    } else if (aSerialElasticLength > (Lsee0 * (1.0 + deltaUseenll))) {
        double deltaFsee0 = getTendonElasticForceRelToFmax() * getMaxIsometricForce();
        double deltaUseel = getTendonElasticLinearStrain();
        Fsee = deltaFsee0 *
               (1.0 + (aSerialElasticLength - Lsee0 * (1.0 + deltaUseenll)) /
                               (deltaUseel * Lsee0));
    }
    return Fsee;
}

double Haeufle2014Muscle::calcIntegralFsee(double aSerialElasticLength) const 
{
    // initialize for case serialElasticLength < lsee,0
    double IntegralFsee = 0.0; 
    double Lsee0 = getTendonSlackLength();
    double deltaUseenll = getTendonElasticNonlinearStrain();
    if (aSerialElasticLength > Lsee0 &&
            aSerialElasticLength <= (Lsee0 * (1.0 + deltaUseenll))) 
    {
        double Fsee = calcFsee(aSerialElasticLength);
        double nuesee = deltaUseenll / getTendonElasticLinearStrain();
        IntegralFsee = Fsee * (aSerialElasticLength - Lsee0) / (nuesee + 1.0);
    } else if (aSerialElasticLength > (Lsee0 * (1.0 + deltaUseenll))) 
    {
        double Fsee = calcFsee(aSerialElasticLength);
        double nuesee = deltaUseenll / getTendonElasticLinearStrain();
        double deltaFsee0 = getTendonElasticForceRelToFmax() * getMaxIsometricForce();
        double deltaUseel = getTendonElasticLinearStrain();
        IntegralFsee = deltaFsee0 * deltaUseenll * Lsee0 / (nuesee + 1.0) +
                       (aSerialElasticLength - Lsee0 * (deltaUseenll + 1.0) / 2.0 *
                                                       (Fsee + deltaFsee0));
    }
    return IntegralFsee;
}

double Haeufle2014Muscle::calcFsde(double serialElasticLengthVelocity,
        double cosPenAng, double parallelElasticForce,
        double contractionForce) const { // contractionForce == F_CE

    double Dse = getTendonMaximumDampingDse();
    double Rse = getTendonOffsetDampingRse();
    double Fmax = getMaxIsometricForce();
    double lceopt = getOptimalFiberLength();
    double Arel0 = getFibreActiveForceVelocityArel0();
    double Brel0 = getFibreActiveForceVelocityBrel0();
    double Fsde = Dse * Fmax * Arel0 / (lceopt * Brel0) *
            ((1.0 - Rse) * cosPenAng *
                            (parallelElasticForce + contractionForce) / Fmax +
                    Rse) *
                  serialElasticLengthVelocity;
    return Fsde;
}

double Haeufle2014Muscle::calcFpde(
    double lengthVelocity, double parallelElasticForce) const 
{
    double Dpe = getFibreMaximumDampingDpe();
    double Rpe = getFibreOffsetDampingRpe();
    double Fmax = getMaxIsometricForce();
    double lceopt = getOptimalFiberLength();
    double Arel0 = getFibreActiveForceVelocityArel0();
    double Brel0 = getFibreActiveForceVelocityBrel0();
    double Fpde = Dpe * Fmax * Arel0 / (lceopt * Brel0) *
                  ((1.0 - Rpe) * parallelElasticForce / Fmax + Rpe) * lengthVelocity;

    return Fpde;
}

double Haeufle2014Muscle::calcC2dash(
        double cosPenAngle, double activation,
        double Fpee, double Brel, double Arel) const 
{
    double Arel0 = getFibreActiveForceVelocityArel0();
    double Brel0 = getFibreActiveForceVelocityBrel0();
    double Dse = getTendonMaximumDampingDse();
    double Rse = getTendonOffsetDampingRse();
    double Fmax = getMaxIsometricForce();
    double Dpe = getFibreMaximumDampingDpe();
    double Rpe = getFibreOffsetDampingRpe();
    double C2woFmax =
            -Dse * Arel0 / (Brel * Brel0 ) *
                    ((1 - Rse) * (Arel - Fpee / Fmax) - Rse/cosPenAngle) +
            cosPenAngle * Dpe * Arel0 / (Brel * Brel0) *
                    ((1 - Rpe) * Fpee / Fmax + Rpe);
    return C2woFmax;
}

double Haeufle2014Muscle::calcC1dash(
        double ldotMTC, double cosPenAngle, double activation, double Fisom, 
        double Fpee, double Fsee, double Arel, double Brel) const 
{
    double Arel0 = getFibreActiveForceVelocityArel0();
    double Brel0 = getFibreActiveForceVelocityBrel0();
    double lceopt = getOptimalFiberLength();
    double Dse = getTendonMaximumDampingDse();
    double Rse = getTendonOffsetDampingRse();
    double Fmax = getMaxIsometricForce();
    double Dpe = getFibreMaximumDampingDpe();
    double Rpe = getFibreOffsetDampingRpe();
    double C1woFmax =
            Dse * Arel0 / Brel0 *
                    (ldotMTC * cosPenAngle / (Brel * lceopt) * (1 - Rse) *
                                    (Arel - Fpee / Fmax) -
                            ldotMTC * Rse / (Brel * lceopt) -
                            (1 - Rse) * (activation * Fisom + Fpee / Fmax) -
                            Rse / cosPenAngle) -
            cosPenAngle * Dpe * Arel0 / Brel0 *
                    ((1 - Rpe) * Fpee / Fmax + Rpe) -
            cosPenAngle * Arel / Brel - Fsee / (Brel * Fmax) +
            cosPenAngle * Fpee / (Brel * Fmax);
    return C1woFmax;
}

double Haeufle2014Muscle::calcC0dash(double ldotMTC,
        double cosPenAngle, double activation, double Fisom, double Fpee,
        double Fsee) const {
    double Arel0 = getFibreActiveForceVelocityArel0();
    double Brel0 = getFibreActiveForceVelocityBrel0();
    double lceopt = getOptimalFiberLength();
    double Dse = getTendonMaximumDampingDse();
    double Rse = getTendonOffsetDampingRse();
    double Fmax = getMaxIsometricForce();
    double C0woFmax =
            Dse * Arel0 * ldotMTC / (lceopt * Brel0) *
                    ((1 - Rse) * cosPenAngle *
                                    (activation * Fisom + Fpee / Fmax) +
                            Rse) -
            cosPenAngle * activation * Fisom - cosPenAngle * Fpee / Fmax +
            Fsee / Fmax;
    return C0woFmax;
}


std::pair<Haeufle2014Muscle::StatusFromInitMuscleState, 
    Haeufle2014Muscle::ValuesFromInitMuscleState>
Haeufle2014Muscle::initMuscleState(const double pathLength,
        const double excitation, const double aSolTolerance,
        const int aMaxIterations) const {
    /** implementing a simple Bisection method to solve the inital fiber length
     * problem:
     * cos(alpha) * ( a * Fisom * Fmax + Fpee) - Fsee = 0
     * within the interval (0, lMTC-lsee0)
     */

     // get interval borders:
    double lower_border = 0.0;
    double upper_border = pathLength - getTendonSlackLength();

    // initialize iteration variable
    int iter = 0;
    double ferr = SimTK::MostPositiveReal; // Solution error

    while ((abs(ferr) > aSolTolerance) && (iter < aMaxIterations))
    {
        double middle_border = (lower_border + upper_border) / 2.0;

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

        double activation_lower = getActivationModel().clampActivation(
            getActivationModel().calculateActivation(excitation, lower_border));
        double activation_middle = getActivationModel().clampActivation(
            getActivationModel().calculateActivation(excitation, middle_border));

        // calculate the total lower and middle values:
        double initalEquilibrium_lower =
                cosPenAng_lower *
                        (activation_lower * Fisom_lower * Fmax + Fpee_lower) -
                Fsee_lower;
        double initalEquilibrium_middle =
                cosPenAng_middle * (activation_middle * Fisom_middle * Fmax +
                                           Fpee_middle) -
                Fsee_middle;
        if ((initalEquilibrium_lower * initalEquilibrium_middle) < 0) {
            upper_border = middle_border;
        } else {
            lower_border = middle_border;
        }
        ferr = upper_border - lower_border;
        iter++;
    }

    // calculate final values:
    double lce_init = ((upper_border + lower_border) / 2.0);
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
