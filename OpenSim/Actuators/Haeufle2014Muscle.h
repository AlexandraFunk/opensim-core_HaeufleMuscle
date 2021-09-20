#ifndef OPENSIM_HAEUFLE2014MUSCLE_H_
#define OPENSIM_HAEUFLE2014MUSCLE_H_
/* -------------------------------------------------------------------------- *
 *                  OpenSim:  Haeufle2014Muscle.h                             *
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
#include <OpenSim/Actuators/osimActuatorsDLL.h>
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
#include <OpenSim/Actuators/MuscleFixedWidthPennationModel.h>

// like the activation dynamic of the fixedwidth pennation models are
#include <OpenSim/Actuators/RockenfellerFirstOrderActivationDynamicModel.h>

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
Haufle et al.\ (2013). Abstract

Hill-type muscle models are commonly used in biomechanical simulations to
predict passive and active muscle forces. Here, a model is presented which
consists of four elements: a contractile element with force-length and
force-velocity relations for concentric and eccentric contractions, a parallel
elastic element, a series elastic element, and a serial damping element. With
this, it combines previously published effects relevant for muscular
contraction, i.e. serial damping and eccentric force-velocity relation. The
model is exemplarily applied to arm movements. The more realistic representation
of the eccentric force-velocity relation results in human-like elbow-joint
flexion.


<B>Reference</B>

D.F.B. Haeufle, M. Guenther, A. Bayer, S. Schmitt (2014) Hill-type muscle model 
with serial damping and eccentric force–velocity relation. Journal of
Biomechanics https://doi.org/10.1016/j.jbiomech.2014.02.009.

@author Maria Hammer
@author Mike Spahr
*/
class OSIMACTUATORS_API Haeufle2014Muscle : public Muscle {
    OpenSim_DECLARE_CONCRETE_OBJECT(Haeufle2014Muscle, Muscle);
public:
//==============================================================================
// PROPERTIES
//==============================================================================
    OpenSim_DECLARE_PROPERTY(default_calcium_concentration, double,
            "Assumed initial activation level if none is assigned.");
    OpenSim_DECLARE_PROPERTY(default_fiber_length, double,
            "Assumed initial fiber length if none is assigned.");

    OpenSim_DECLARE_PROPERTY(fibre_exponent_descending_active_force_length, double,
            "The exponent of the normalized bell curve in its descending"
            "limb as used in Haeufle et al (2014).");
    OpenSim_DECLARE_PROPERTY(fibre_width_descending_active_force_length, double,
            "The width of the normalized bell curve in its descending"
            "limb as used in Haeufle et al (2014).");
    OpenSim_DECLARE_PROPERTY(fibre_exponent_ascending_active_force_length, double,
            "The exponent of the normalized bell curve in its ascending"
            "limb as used in Haeufle et al (2014).");
    OpenSim_DECLARE_PROPERTY(fibre_width_ascending_active_force_length, double,
            "The width of the normalized bell curve in its ascending"
            "limb as used in Haeufle et al (2014).");
    OpenSim_DECLARE_PROPERTY(fibre_active_force_velocity_arel0, double,
            "Derived from the classical Hill constant a");
    OpenSim_DECLARE_PROPERTY(fibre_active_force_velocity_brel0, double,
            "Derived from the classical Hill constant b");
    OpenSim_DECLARE_PROPERTY(fibre_max_eccentric_force_rel_to_fmax, double,
            "Maximal Force during eccentric lengthening");
    OpenSim_DECLARE_PROPERTY(fibre_eccentric_slopefactor, double,
            "The ratio between the eccentric and concentric derivatives "
            "dF/dVce");
    OpenSim_DECLARE_PROPERTY(fibre_elastic_zero_length, double,
            "Zero length of the parallel elastic element");
    OpenSim_DECLARE_PROPERTY(fibre_elastic_exponent, double,
            "Exponent of the parallel elastic element");
    OpenSim_DECLARE_PROPERTY(fibre_elastic_force_rel_to_fmax, double,
            " Parallel elastic element force relative to Fmax at "
            "Lceopt*(1+dWdes)");
    OpenSim_DECLARE_PROPERTY(
            tendon_elastic_nonlinear_strain, double, "");
    OpenSim_DECLARE_PROPERTY(tendon_elastic_linear_strain, double, "");
    OpenSim_DECLARE_PROPERTY(tendon_elastic_force_rel_to_fmax, double, "");
    OpenSim_DECLARE_PROPERTY(tendon_maximum_damping_dse, double,
            "The dse damping factor as used in Moerl et al (2012)");
    OpenSim_DECLARE_PROPERTY(tendon_offset_damping_rse, double,
            "The rse damping factor as used in Moerl et al (2012)");
    OpenSim_DECLARE_PROPERTY(fibre_maximum_damping_dpe, double,
            "The dpe damping factor as used in Moerl et al (2012)");
    OpenSim_DECLARE_PROPERTY(fibre_offset_damping_rpe, double,
            "The rpe damping factor as used in Moerl et al (2012)");
    OpenSim_DECLARE_PROPERTY(maximum_pennation_angle, double,
            "Maximum pennation angle (in radians).");
    OpenSim_DECLARE_PROPERTY(activation_Hatze_time_constant, double,
            "Time constant, in 1/seconds. (overridden when this is a "
            "subcomponent of a Muscle)");
    OpenSim_DECLARE_PROPERTY(activation_exponent, double,
            "Hatze Coefficient (overridden when this is a subcomponent of a "
            "Muscle)");
    OpenSim_DECLARE_PROPERTY(activation_optimal_calcium_concentration_fraction,
            double, "roh_0 * gamma_C");
    OpenSim_DECLARE_PROPERTY(activation_minimum, double,
            "Activation lower bound (Hatze constant 0.005) equal to Kuh0");

//==============================================================================
// OUTPUTS
//==============================================================================
    OpenSim_DECLARE_OUTPUT(
            Fpee, double, getFpee, SimTK::Stage::Dynamics);
    OpenSim_DECLARE_OUTPUT(
            Fpde, double, getFpde, SimTK::Stage::Dynamics);
    OpenSim_DECLARE_OUTPUT(
            Fsee, double, getFsee, SimTK::Stage::Dynamics);
    OpenSim_DECLARE_OUTPUT(
            Fsde, double, getFsde, SimTK::Stage::Dynamics);

//==============================================================================
// CONSTRUCTORS
//==============================================================================
    /** Default constructor. Produces a non-functional empty muscle. */
    Haeufle2014Muscle();

     /** Constructs a functional muscle using default curves and activation model
    parameters. The tendon is assumed to be elastic, full fiber dynamics are
    solved, and activation dynamics are included.
        @param aName The name of the muscle.
        @param aMaxIsometricForce The force generated by the muscle when fully
    activated at its optimal resting length with a contraction velocity of zero.
        @param aOptimalFiberLength The optimal length of the muscle fiber.
        @param aTendonSlackLength The resting length of the tendon.
        @param aPennationAngle The angle of the fiber (in radians) relative to
    the tendon when the fiber is at its optimal resting length. */
    Haeufle2014Muscle(const std::string& aName, double aMaxIsometricForce,
            double aOptimalFiberLength, double aTendonSlackLength,
            double aPennationAngle);


//==============================================================================
// GET METHODS
//==============================================================================

    /** @returns The default activation level that is used as an initial
    condition if none is provided by the user. */
    double getDefaultCalciumConcentration() const;

    /** @returns The default fiber length that is used as an initial condition
    if none is provided by the user. */
    double getDefaultFiberLength() const;

    /** @returns The exponent of the descending limb of the normalized
    bell curve */
    double getFibreExponentDescendingActiveForceLength() const;

    /** @returns The width of the descending limb of the normalized
    bell curve */
    double getFibreWidthDescendingActiveForceLength() const;

    /** @returns The exponent of the ascending limb of the normalized
    bell curve */
    double getFibreExponentAscendingActiveForceLength() const;

    /** @returns The width of the ascending limb of the normalized
    bell curve */
    double getFibreWidthAscendingActiveForceLength() const;

    /** @returns The concentric contraction parameter which is called Arel in
    Haeufle et al */
    double getFibreActiveForceVelocityArel0() const;

    /** @returns The concentric contraction parameter which is called Brel in
    Haeufle et al */
    double getFibreActiveForceVelocityBrel0() const;

    /** @returns The maximum force at eccentric lengthening which is called Fe
    in Haeufle et al */
    double getFibreMaxEccentricForceRelToFmax() const;

    /** @returns The fibreEccentricSlopefactor of the curve which is called Se 
     * in Haeufle et al */
    double getFibreEccentricSlopefactor() const;

    /** @returns The parallel elastic elements zero length (L_PEE,0) */
    double getFibreElasticZeroLength() const;

    /** @returns The parallel elastic elements exponent (nue_PEE) */
    double getFibreElasticExponent() const;

    /** @returns The parallel elastic elemnts force relative to Fmax at
        Lceopt*(1+dWdes) */
    double getFibreElasticForceRelToFmax() const;

    /** @returns The relative stretch at nonlinear/linear transition in Fsee */
    double getTendonElasticNonlinearStrain() const;

    /** @returns The relative stretch in the linear part for force increase
     * delta Fsee,0 */
    double getTendonElasticLinearStrain() const;

    /** @returns The force at nonlinear/linear transition in Fsee */
    double getTendonElasticForceRelToFmax() const;

    /** @returns The dse damping factor */
    double getTendonMaximumDampingDse() const;

    /** @returns The rse damping factor */
    double getTendonOffsetDampingRse() const;

    /** @returns The dse damping factor */
    double getFibreMaximumDampingDpe() const;

    /** @returns The rse damping factor */
    double getFibreOffsetDampingRpe() const;

    /** @returns The MuscleFixedWidthPennationModel owned by this model. */
    const MuscleFixedWidthPennationModel& getPennationModel() const;

    /** @returns The RockenfellerFirstOrderActivationDynamicModel owned by this
    model. */
    const RockenfellerFirstOrderActivationDynamicModel& 
    getActivationModel() const;

    /** @returns The current calcium concentration derivative value */
    double getCalciumConcentrationDerivative(const SimTK::State& s) const;

    /** @returns The current calcium concentration */
    double getCalciumConcentration(const SimTK::State& s) const;

    /** @returns The minimum fiber length **/
    double getMinimumFiberLength() const;

    /** @returns The minimum fiber length along the tendon **/
    double getMinimumFiberLengthAlongTendon() const;

    /** @returns The parallel elastic elements force **/
    double getFpee(const SimTK::State& s) const;

    /** @returns The parallel damping elements force **/
    double getFpde(const SimTK::State& s) const;

    /** @returns The serial elastic elements force **/
    double getFsee(const SimTK::State& s) const;

    /** @returns The serial damping elements force **/
    double getFsde(const SimTK::State& s) const;

//==============================================================================
// SET METHODS
//==============================================================================

    /** @param calciumConcentration The default concentration level that is used 
     * to initialize the muscle. */
    void setDefaultCalciumConcentration(double calciumConcentration);

    /** @param fiberLength The default fiber length that is used to initialize
    the muscle. */
    void setDefaultFiberLength(double fiberLength);

    /**
    @param fibreExponentDescendingActiveForceLength
        The exponent of the descending limb of the normalized
        bell curve.
    @param fibreWidthDescendingActiveForceLength
        The width of the descending limb of the normalized
        bell curve.
    @param fibreExponentAscendingActiveForceLength
        The exponent of the ascending limb of the normalized
        bell curve.
    @param fibreWidthAscendingActiveForceLength
        The width of the ascending limb of the normalized
        bell curve.
    */
    void setActiveForceLengthParameters(
            double fibreExponentDescendingActiveForceLength,
            double fibreWidthDescendingActiveForceLength,
            double fibreExponentAscendingActiveForceLength,
            double fibreWidthAscendingActiveForceLength);

    /** 
    @param aFibreActiveForceVelocityArel0
        The parameter for the concentric contraction dynamics of the CE Arel0
    @param aFibreActiveForceVelocityBrel0
        The parameter for the concentric contraction dynamics of the CE Brel0
    @param aFibreMaxEccentricForceRelToFmax
        The force value at a maximal eccentric extension as defined by van 
        Soest & Bobbert, 1993
    @param aFibreEccentricSlopefactor
        The corresponding slope factor as defined by van Soest & Bobbert, 1993
    */
    void setFiberVelocityParameters(double aFibreActiveForceVelocityArel0,
            double aFibreActiveForceVelocityBrel0,
            double aFibreMaxEccentricForceRelToFmax, double aFibreEccentricSlopefactor);

    /**
    @param aFibreElasticZeroLength
        The parallel elastic elements zero lenght
    @param aFibreElasticExponent
        The parallel elastic elements exponent
    @param aFibreElasticForceRelToFmax
        The parallel elastic elemnsts force relative to Fmax at Lceopt*(1+dWdes)
    */
    void setParallelElasticParameters(double aFibreElasticZeroLength,
            double aFibreElasticExponent,
            double aFibreElasticForceRelToFmax);

    /**
    @param aTendonElasticForceRelToFmax
        The force at non-linear/linear transition which is muscle specific
    */
    void setTendonElasticForceRelToFmax(
            double aTendonElasticForceRelToFmax);

    /**
    @param aForceAtNonlinearLinearTransition
        The force at non-linear/linear transition which is muscle specific
    */
    void setTendonElasticNonlinearStrain(
            double aTendonElasticNonlinearStrain);

    /**
    @param aForceAtNonlinearLinearTransition
        The force at non-linear/linear transition which is muscle specific
    */
    void setTendonElasticLinearStrain(double aTendonElasticLinearStrain);

    /**
    @param aTendonMaximumDampingDse
        The Dse damping factor
    @param aTendonOffsetDampingRse
        The Rse damping factor
    */
    void setTendonDampingParams(double aTendonMaximumDampingDse,
            double aTendonOffsetDampingRse);

    /**
    @param aTendonMaximumDampingDse
        The Dse damping factor
    @param aTendonOffsetDampingRse
        The Rse damping factor
    */
    void setParallelDampingParams(
            double aFibreMaximumDampingDpe, double aFibreOffsetDampingRpe);

    /** @param[out] s The state of the system.
        @param fiberLength The desired fiber length (m). */
    void setFiberLength(SimTK::State& s, double fiberLength) const;


//==============================================================================
// MUSCLE.H INTERFACE
//==============================================================================
    /**  This is neccessary  since it is implemented in the muscle.h interface
     * as virtual but sets here the normalized calcium concentration and not the
     * activation */
    virtual void setActivation(SimTK::State& s, double activation) const;

    // End of Muscle's State Dependent Accessors.
    //@}

    /** Actuator interface for a muscle computes the tension in the muscle
        and applied by the tendon to bones (i.e. not the fiber force) */
    double computeActuation(const SimTK::State& s) const override;

    /** Computes the fiber length such that the fiber and tendon are developing
    the same force, distributing the velocity of the entire musculotendon
    actuator between the fiber and tendon according to their relative
    stiffnesses. Assuming muscle-tendon velocity as zero
    @param[in,out] s The state of the system.
    @throws MuscleCannotEquilibrate
    */
    void computeInitialFiberEquilibrium(SimTK::State& s) const override;

    /** Adjust the properties of the muscle after the model has been scaled. The
    optimal fiber length and tendon slack length are each multiplied by the
    ratio of the current path length and the path length before scaling. */
    void extendPostScale(
            const SimTK::State& s, const ScaleSet& scaleSet) override;

    /*  @returns the fiber force consisting of Fce and Fpee generated by a rigid tendon equilibrium
    muscle model
    @param s the state of the system
    @param aActivation the activation of the muscle */
    double calcInextensibleTendonActiveFiberForce(
            SimTK::State& s, double aActivation) const override final;

protected:

    /** Gets the derivative of an actuator state by index.
    @param s The state.
    @param aStateName The name of the state to get.
    @return The value of the state derivative. */
    double getStateVariableDeriv(
            const SimTK::State& s, const std::string& aStateName) const;

    /** Sets the derivative of an actuator state specified by name.
        @param s The state.
        @param aStateName The name of the state to set.
        @param aValue The value to which the state should be set. */
    void setStateVariableDeriv(const SimTK::State& s,
            const std::string& aStateName, double aValue) const;


    //--------------------------------------------------------------------------
    // HAEUFLE MUSCLE CALCULATION FUNCTIONS FOR ALL NECESSARY FORCES
    //--------------------------------------------------------------------------

    /** Evaluates the active-force-length curve at a fiber
    length of 'fiberLength'. */
    double calcFisom(double fiberLength) const;

    /** Evaluates the force-velocity curve at a fiber velocity of
    'fiberVelocity'. */
    double calcNormFce(double fiberVelocity,
            double fiberLength, double activation) const;

    /** Calculate the fiber length and activation dependent normalized Hill
    parameter Arel for the conccentric case */
    double calcArel(double fiberLength, double activation, double Fisom) const;

    /** Calculate the activation dependent normalized Hill
    parameter Brel for the conccentric case */
    double calcBrel(double activation) const;

    /** Calculate the activation and Fisom dependent
    normalized Hill parameter Arele for the eccentric case */
    double calcArele(double activation, double Fisom) const;

    /** Calculate the activation and correspondingForceLengthValue dependent
    normalized Hill parameter Brele for the eccentric case */
    double calcBrele(double fiberLength, double activation,
            double Fisom) const;

    /** Calculates the parallel elastic force at a fiber length of
     * 'FiberLength'. */
    double calcFpee(double fiberLength) const;

    /** Calculates the integral of the parallel elastic force at 
     * a fiber length of 'FiberLength'. */
    double calcIntegralFpee(double fiberLength) const;

    /** Calculates the tendon-force-length force at a tendon length of
    'aSerialElasticLength'. */
    double calcFsee(double aSerialElasticLength) const;

    /** Calculates the integral of the tendon-force-length force at a 
     * tendon length of 'aSerialElasticLength'. */
    double calcIntegralFsee(double aSerialElasticLength) const;

    /** Calculates the tendon damping force Fsde which depends on the serial
     * elastic length velocity and the elasticTendonForce (Fsee) */
    double calcFsde(double serialElasticLengthVelocity, double cosPenAng,
            double parallelElasticForce,
            double contractionForce) const;
 
    /** Calculates the parallel damping force Fpde which depends on the
     *  length velocity and the parallelElasticForce (Fpee) */
    double calcFpde(double lengthVelocity, double parallelElasticForce) const;

    /** Helper function to calculate C2 dashed which is needed during the 
     * calculation of the fiber velocity 
     */
    double calcC2dash(double cosPenAngle, double activation, double Fpee,
            double Brel, double Arel) const;

    /** Helper function to calculate C1 dashed which is needed during the 
     * calculation of the fiber velocity 
     */
    double calcC1dash(double ldotMTC, double cosPenAngle, double activation,
            double Fisom, double Fpee, double Fsee, double Arel,
            double Brel) const;

    /** Helper function to calculate C0 dashed which is needed during the 
     * calculation of the fiber velocity 
     */
    double calcC0dash(double ldotMTC, double cosPenAngle,
            double activation, double Fisom, double Fpee, double Fsee) const;

    /** Calculate the maximum contraction velocity at the beginning of the 
     * simulation.*/
    double calcVmax() const;

    //--------------------------------------------------------------------------
    // CALCULATIONS
    //--------------------------------------------------------------------------
    /** @name Muscle State Dependent Calculations
     *  Developers must override these methods to implement the desired behavior
     *  of their muscle models. Unless you are augmenting the behavior
     *  of an existing muscle class or writing a new derived class, you do not
     *  have access to these methods.
     */
    //@{
    /** calculate muscle's position related values such fiber and tendon
       lengths, normalized lengths, pennation angle, etc... */
    void calcMuscleLengthInfo(
            const SimTK::State& s, MuscleLengthInfo& mli) const override;

    /** calculate muscle's fiber velocity and pennation angular velocity, etc...
     */
    void calcFiberVelocityInfo(
            const SimTK::State& s, FiberVelocityInfo& fvi) const override;

    /** calculate muscle's active and passive force-length, force-velocity,
        tendon force, relationships and their related values */
    void calcMuscleDynamicsInfo(
            const SimTK::State& s, MuscleDynamicsInfo& mdi) const override;

    /** calculate muscle's fiber and tendon potential energy */
    void calcMusclePotentialEnergyInfo(
            const SimTK::State& s, MusclePotentialEnergyInfo& mpei) const override;

//==============================================================================
// MODELCOMPONENT INTERFACE REQUIREMENTS
//==============================================================================
    /** Sets up the ModelComponent from the model, if necessary */
    void extendConnectToModel(Model& model) override;

    /** Creates the ModelComponent so that it can be used in simulation */
    void extendAddToSystem(SimTK::MultibodySystem& system) const override;

    /** Initializes the state of the ModelComponent */
    void extendInitStateFromProperties(SimTK::State& s) const override;

    /** Sets the default state for the ModelComponent */
    void extendSetPropertiesFromState(const SimTK::State& s) override;

    /** Computes and sets state variable derivatives */
    void computeStateVariableDerivatives(const SimTK::State& s) const override;

private:
    // The name used to access the calcium concentration state.
    static const std::string STATE_CALCIUM_CONCENTRATION;
    // The name used to access the fiber length state.
    static const std::string STATE_FIBER_LENGTH_NAME;

    void setNull();
    void constructProperties();

    // Rebuilds muscle model if any of its properties have changed.
    void extendFinalizeFromProperties() override;

    // Subcomponents owned by the muscle. The properties of these subcomponents
    // are set in extendFinalizeFromProperties() from the properties of the
    // muscle.
    MemberSubcomponentIndex penMdlIdx{
            constructSubcomponent<MuscleFixedWidthPennationModel>("penMdl")};
    MemberSubcomponentIndex actMdlIdx{
            constructSubcomponent<RockenfellerFirstOrderActivationDynamicModel>(
                    "actMdl")};

    // Status flag returned by ().
    enum StatusFromInitMuscleState {
        Success_Converged,
        Failure_MaxIterationsReached
    };

    // Associative array of values returned by initMuscleState():
    // solution_error, iterations, fiber_length, fiber_velocity,
    // tendon_force and activation
    typedef std::map<std::string, double> ValuesFromInitMuscleState;


    // Returns the maximum of the minimum fiber length and the current fiber
    // length.
    double clampFiberLength(double lce) const;

    // The tendon Damping Force which is not included in the Muscle.h 
    // MuscleDynamicsInfo struct but is necessary for the Haeufle2014Muscle 
    // modell.
    double m_tendonDampingForce;

    double m_minimumFiberLength;
    double m_minimumFiberLengthAlongTendon;

    std::pair<StatusFromInitMuscleState, ValuesFromInitMuscleState>
    initMuscleState(const double pathLength, const double excitation,
            const double aSolTolerance, const int aMaxIterations) const;

};

} //end of namespace OpenSim

#endif // OPENSIM_HAEUFLE2014MUSCLE_H_
