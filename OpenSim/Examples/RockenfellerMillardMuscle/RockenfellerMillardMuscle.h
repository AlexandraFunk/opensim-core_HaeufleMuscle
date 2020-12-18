#ifndef OPENSIM_ROCKENFELLER_MILLARD_MUSCLE_H_
#define OPENSIM_ROCKENFELLER_MILLARD_MUSCLE_H_
/* -------------------------------------------------------------------------- *
 *                        OpenSim:  RockenfellerMillardMuscle.h                         *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Copyright (c) 2005-2017 Stanford University and the Authors                *
 * Author(s): Peter Loan, Ajay Seth                                           *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied    *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */


// INCLUDE
#include <OpenSim/OpenSim.h>

namespace OpenSim {

//=============================================================================
//=============================================================================
/**
 * This class extends a Millard2012EquilibriumMuscle by including three
 * additional states to model the fatigue and recovery of muscle fibers.
 * The equations for these states are (loosely) based on the following paper:
 * Liu, Jing Z., Brown, Robert W., Yue, Guang H., "A Dynamical Model of Muscle
 * Activation, Fatigue, and Recovery," Biophysical Journal, Vol. 82, Issue 5,
 * pp. 2344-2359, 2002.
 *
 * NOTE: The primary purpose of this muscle model is to serve as an example
 *       for extending existing OpenSim muscle models and it is not intended
 *       for research. The implementation of the fatigue dynamics have not
 *       been tested.
 *
 * @author Ajay Seth (based on Millard2012EquilibriumMuscle)
 * @contributor Peter Loan (originally based on Thelen2003Muscle)
 *
 * The Muscle base class, specifies the interface that must be implemented
 * by derived muscle classes. The FatigableMuscle derives from
 * Millard2012EquilibriumMuscle, which is a concrete implementation of the
 * Muscle interface. The dynamics for fatigue are added by overriding methods
 * extendAddToSystem() which allocates the additional states and
 * computeStateVariableDerivatives() to specify their dynamics (derivatives).
 *
 * @see Millard2012EquilibriumMuscle
 * @see Muscle
 */
class RockenfellerMillardMuscle : public Millard2012EquilibriumMuscle {
OpenSim_DECLARE_CONCRETE_OBJECT(RockenfellerMillardMuscle,
                                Millard2012EquilibriumMuscle);
public:
//=============================================================================
// PROPERTIES
//=============================================================================
    // OpenSim_DECLARE_PROPERTY(fatigue_factor, double,
    //     "percentage of active motor units that fatigue in unit time");
    // OpenSim_DECLARE_PROPERTY(recovery_factor, double,
    //     "percentage of fatigued motor units that recover in unit time");
    // OpenSim_DECLARE_PROPERTY(default_target_activation, double,
    //     "default state value used to initialize the muscle's target activation");
    // OpenSim_DECLARE_PROPERTY(default_active_motor_units, double,
    //     "default state value for the fraction of motor units that are active");
    // OpenSim_DECLARE_PROPERTY(default_fatigued_motor_units, double,
    //     "default state value for the fraction of motor units that are fatigued");

    OpenSim_DECLARE_PROPERTY(time_constant_hatze, double,
        "Time constant, in 1/seconds.");
    OpenSim_DECLARE_PROPERTY(nue, double,
        "Hatze Coefficient");
    OpenSim_DECLARE_PROPERTY(roh_0, double,
        "Hatze constant [l/mol] from Rockenfeller2018");
    OpenSim_DECLARE_PROPERTY(gamma_C, double,
        "Hatze Coefficient"); // TODO add description
    // OpenSim_DECLARE_PROPERTY(default_normalized_calcium_concentration, double,
    //     "Normalized Calcium concentration (by Ca2+_max).");

public:
//=============================================================================
// METHODS
//=============================================================================
    //-------------------------------------------------------------------------
    // CONSTRUCTION
    //-------------------------------------------------------------------------
    RockenfellerMillardMuscle();
    RockenfellerMillardMuscle(const std::string &name, double maxIsometricForce,
                       double optimalFiberLength, double tendonSlackLength,
                       double pennationAngle, double time_constant_hatze,
                       double nue, double roh_0, double gamma_C);

    // employs the default destructor, copy constructor and copy assignment
    // that are automatically supplied by the compiler if none are defined

    //-------------------------------------------------------------------------
    // GET & SET Properties
    //-------------------------------------------------------------------------
    double getTimeConstantHatze() const { return get_time_constant_hatze(); }
    void setTimeConstantHatze(double aTimeConstantHatze);
    double getNue() const { return get_nue(); }
    void setNue(double aNue);
    double getRoh_0() const { return get_roh_0(); }
    void setRoh_0(double roh_0);
    double getGamma_C() const { return get_gamma_C(); }
    void setGamma_C(double gamma_C);

    /** default values for states */
    // double getDefaultNormalizedCalciumConcentration() const {
    //     return get_default_normalized_calcium_concentration();
    // }
    // void setDefaultNormalizedCalciumConcentration(double anormalized_calcium_concentration);

    //-------------------------------------------------------------------------
    // GET & SET States and their derivatives
    //-------------------------------------------------------------------------
    // double getNormalizedCalciumConcentration(const SimTK::State& s) const;
    // void setNormalizedCalciumConcentration(SimTK::State& s, double normalizedCa) const;

    // double getNormalizedCalciumConcentrationDeriv(const SimTK::State& s) const;
    // void setNormalizedCalciumConcentrationDeriv(const SimTK::State& s,
    //                                             double normalizedCaDeriv) const;

    double rho(double lce) const;

protected:
    // Model Component Interface
    /** add new dynamical states to the multibody system corresponding
        to this muscle */
    void extendAddToSystem(SimTK::MultibodySystem& system) const override;
    /** initialize muscle state variables from properties. For example, any
        properties that contain default state values */
    void extendInitStateFromProperties(SimTK::State& s) const override;
    /** use the current values in the state to update any properties such as
        default values for state variables */
    void extendSetPropertiesFromState(const SimTK::State& s) override;

    //-------------------------------------------------------------------------
    // COMPUTATIONS
    //-------------------------------------------------------------------------
    /** Compute the derivatives for state variables added by this muscle */
    void computeStateVariableDerivatives(const SimTK::State& s) const override;

    /** Calculate the dynamics-related values associated with the muscle state
    (from the active- and passive-force-length curves, the force-velocity curve,
    and the tendon-force-length curve). The last entry is a SimTK::Vector
    containing the passive conservative (elastic) fiber force and the passive
    non-conservative (damping) fiber force. */
    void calcMuscleDynamicsInfo(const SimTK::State& s,
                                MuscleDynamicsInfo& mdi) const override;

private:
    /** construct the new properties and set their default values */
    void constructProperties();

    /*  @param fiso the maximum isometric force the fiber can generate
    @param a activation
    @param fal the fiber active-force-length multiplier
    @param fv the fiber force-velocity multiplier
    @param fpe the fiber force-length multiplier
    @param dlceN the normalized fiber velocity
    @param cosphi the cosine of the pennation angle
    @returns Vec4
        [0] total fiber force
        [1] active fiber force
        [2] conservative passive fiber force
        [3] non-conservative passive fiber force */
    SimTK::Vec4 calcFiberForce(double fiso, double a, double fal, double fv,
            double fpe, double dlceN) const;

    /*  @param fiso the maximum isometric force the fiber can generate
    @param a activation
    @param fal the fiber active-force-length multiplier
    @param fv the fiber force-velocity multiplier
    @param fpe the fiber force-length multiplier
    @param sinphi the sine of the pennation angle
    @param cosphi the cosine of the pennation angle
    @param lce the fiber length
    @param lceN the normalized fiber length
    @param optFibLen the optimal fiber length
    @returns the stiffness of the fiber in the direction of the fiber */
    double calcFiberStiffness(double fiso, double a, double fv, double lceN,
            double optFibLen) const;


    /*  @param fiberForce the force, in Newtons, developed by the fiber
    @param fiberStiffness the stiffness, in N/m, of the fiber
    @param lce the fiber length
    @param sinphi the sine of the pennation angle
    @param cosphi the cosine of the pennation angle
    @returns the partial derivative of fiber force along the tendon with
    respect to small changes in fiber length (in the direction of the fiber) */
    double calc_DFiberForceAT_DFiberLength(double fiberForce,
            double fiberStiffness, double lce, double sinPhi,
            double cosPhi) const;

    /*  @param dFm_d_lce stiffness of the fiber in the direction of the fiber
    @param sinphi the sine of the pennation angle
    @param cosphi the cosine of the pennation angle
    @param lce the fiber length
    @returns the stiffness of the fiber in the direction of the tendon */
    double calc_DFiberForceAT_DFiberLengthAT(
            double dFm_d_lce, double sinPhi, double cosPhi, double lce) const;

//=============================================================================
};

} // end of namespace OpenSim

#endif
