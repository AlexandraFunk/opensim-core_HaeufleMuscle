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
#include <OpenSim/Actuators/MuscleFixedWidthPennationModel.h>

#include "RockenfellerFirstOrderActivationDynamicModel.h"
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
class OSIMPLUGIN_API Haeufle2014Muscle : public Muscle {
    OpenSim_DECLARE_CONCRETE_OBJECT(Haeufle2014Muscle, Muscle);
public:
//==============================================================================
// PROPERTIES
//==============================================================================
    // TODO check which other properties are necessary?
    // which of these properties are essential to describe a muscle and are needed for calculations?
    OpenSim_DECLARE_PROPERTY(default_activation, double,
            "Assumed initial activation level if none is assigned.");
    OpenSim_DECLARE_PROPERTY(default_fiber_length, double,
            "Assumed initial fiber length if none is assigned.");

    OpenSim_DECLARE_UNNAMED_PROPERTY(
            HaeufleActiveForceLengthCurve, "Active-force-length curve.");
    OpenSim_DECLARE_UNNAMED_PROPERTY(
            HaeufleForceVelocityCurve, "Force-velocity curve.");
    OpenSim_DECLARE_UNNAMED_PROPERTY(
            HaeufleFiberForceLengthCurve, "Passive-force-length curve.");
    OpenSim_DECLARE_UNNAMED_PROPERTY(
            HaeufleTendonForceLengthCurve, "Tendon-force-length curve.");
    OpenSim_DECLARE_UNNAMED_PROPERTY(
            HaeufleTendonDampingCurve, "Tendon-damping curve");


//==============================================================================
// OUTPUTS
//==============================================================================


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
    double getDefaultActivation() const;

    /** @returns The default fiber length that is used as an initial condition
    if none is provided by the user. */
    double getDefaultFiberLength() const;


    /** @returns The ActiveForceLengthCurve used by this model. */
    const HaeufleActiveForceLengthCurve& getActiveForceLengthCurve() const;
    /** @returns The ForceVelocityCurve used by this model. */
    const HaeufleForceVelocityCurve& getForceVelocityCurve() const;
    /** @returns The FiberForceLengthCurve used by this model. */
    const HaeufleFiberForceLengthCurve& getFiberForceLengthCurve() const;
    /** @returns The TendonForceLengthCurve used by this model. */
    const HaeufleTendonForceLengthCurve& getTendonForceLengthCurve() const;
    /** @returns The TendonDampingCurve used by this model. */
    const HaeufleTendonDampingCurve& getTendonDampingCurve() const;

    /** @returns The MuscleFixedWidthPennationModel owned by this model. */
    const MuscleFixedWidthPennationModel& getPennationModel() const;

    /** @returns The RockenfellerFirstOrderActivationDynamicModel owned by this
    model. */
    const RockenfellerFirstOrderActivationDynamicModel&
    getActivationModel() const;



    // TODO add more getter methods which are needed during calculation, for example
    // getFiberVelocity, getActivation (see Muscle.h) ....




//==============================================================================
// SET METHODS
//==============================================================================

    /** @param activation The default activation level that is used to
    initialize the muscle. */
    void setDefaultActivation(double activation);

    /** @param fiberLength The default fiber length that is used to initialize
    the muscle. */
    void setDefaultFiberLength(double fiberLength);


    /** @param aActiveForceLengthCurve The HaeufleActiveForceLengthCurve used by the
    muscle model to scale active fiber force as a function of fiber length. */
    void setActiveForceLengthCurve(
            HaeufleActiveForceLengthCurve& aActiveForceLengthCurve);

    /** @param aForceVelocityCurve The HaeufleForceVelocityCurve used by the muscle
    model to calculate the derivative of fiber length. */
    void setForceVelocityCurve(HaeufleForceVelocityCurve& aForceVelocityCurve);

    /** @param aFiberForceLengthCurve The HaeufleFiberForceLengthCurve used by the
    muscle model to calculate the passive force the muscle fiber generates as a
    function of fiber length. */
    void setFiberForceLengthCurve(
            HaeufleFiberForceLengthCurve& aFiberForceLengthCurve);

    /** @param aTendonForceLengthCurve The HaeufleTendonForceLengthCurve used by the
    muscle model to calculate the force exerted by the tendon as a function of
    tendon length. */
    void setTendonForceLengthCurve(
            HaeufleTendonForceLengthCurve& aTendonForceLengthCurve);

    /** @param aTendonDampingCurve The HaeufleTendonDampingCurve used by
    the muscle model. */
    void setTendonDampingCurve(
            HaeufleTendonDampingCurve& aTendonDampingCurve);

    /** @param[out] s The state of the system.
        @param fiberLength The desired fiber length (m). */
    void setFiberLength(SimTK::State& s, double fiberLength) const;


//==============================================================================
// MUSCLE.H INTERFACE
//==============================================================================
    /** DEPRECATED: only for backward compatibility */
    virtual void setActivation(SimTK::State& s, double activation) const;

    // End of Muscle's State Dependent Accessors.
    //@}

    /** Actuator interface for a muscle computes the tension in the muscle
        and applied by the tendon to bones (i.e. not the fiber force) */
    double computeActuation(const SimTK::State& s) const override;

    /** Computes the fiber length such that the fiber and tendon are developing
    the same force, distributing the velocity of the entire musculotendon
    actuator between the fiber and tendon according to their relative
    stiffnesses.
    @param[in,out] s The state of the system.
    @throws MuscleCannotEquilibrate
    */
    void computeInitialFiberEquilibrium(SimTK::State& s) const override {
        computeFiberEquilibrium(s, false);
    }

    /** Computes the fiber length such that the fiber and tendon are developing
        the same force, either assuming muscle-tendon velocity as provided
        by the state or zero as designated by the useZeroVelocity flag.
        @param[in,out] s         The state of the system.
        @param solveForVelocity  Flag indicating to solve for fiber velocity,
                                 which by default is false (zero fiber-velocity)
        @throws MuscleCannotEquilibrate
    */
    void computeFiberEquilibrium(
            SimTK::State& s, bool solveForVelocity = false) const;


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

//==============================================================================
// MUSCLE INTERFACE REQUIREMENTS
//==============================================================================
    /**Developer Access to intermediate values calculate by the muscle model*/
    const MuscleLengthInfo& getMuscleLengthInfo(
                    const SimTK::State& s) const;
    MuscleLengthInfo& updMuscleLengthInfo(const SimTK::State& s) const;

    const FiberVelocityInfo& getFiberVelocityInfo(const SimTK::State& s) const;
    FiberVelocityInfo& updFiberVelocityInfo(const SimTK::State& s) const;

    const MuscleDynamicsInfo& getMuscleDynamicsInfo(
            const SimTK::State& s) const;
    MuscleDynamicsInfo& updMuscleDynamicsInfo(const SimTK::State& s) const;

    const MusclePotentialEnergyInfo& getMusclePotentialEnergyInfo(
            const SimTK::State& s) const;
    MusclePotentialEnergyInfo& updMusclePotentialEnergyInfo(
            const SimTK::State& s) const;

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

    /** Computes state variable derivatives */
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

};

} //end of namespace OpenSim

#endif // OPENSIM_Haeufle2014Muscle_h__
