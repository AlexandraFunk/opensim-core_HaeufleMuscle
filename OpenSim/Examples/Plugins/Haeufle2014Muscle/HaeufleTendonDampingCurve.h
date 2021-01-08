#ifndef OPENSIM_HAEUFLE_TENDON_DAMPING_CURVE_H_
#define OPENSIM_HAEUFLE_TENDON_DAMPING_CURVE_H_

/* -------------------------------------------------------------------------- *
 *                     OpenSim:  HaeufleTendonDampingCurve.h                  *
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
#    include <OpenSim/Actuators/osimActuatorsDLL.h>
#    include <OpenSim/Simulation/Model/ModelComponent.h>

#ifdef SWIG
	#ifdef OSIMACTUATORS_API
		#undef OSIMACTUATORS_API
		#define OSIMACTUATORS_API
	#endif
#endif

namespace OpenSim {
/** This class serves as a TendonDampingCurve as described by Haeufle et
    al in Hill-type muscle model with serial damping and eccentric force–velocity
    relation. Journal of Biomechanics (D.F.B. Haeufle, M. Guenther, A. Bayer, S.
    Schmitt (2014))
    This class depends on the following parameters:

    @param max_isometric_force
    @param optimal_fiber_length
    @param dse_damping_factor
    @param rse_damping_factor
    @param concentric_contraction_a_rel0
    @param concentric_contraction_b_rel0


    <B>Default Parameter Values</B>
    \verbatim
    maxIsometricForce ............ 1500
    optimalFiberLength ........... 0.12
    dseDampingFactor ............. 0.3
    rseDampingFactor ............. 0.01
    concentricContractionARel0 ... 0.2
    concentricContractionBRel0 ... 2.0
    \endverbatim

    <B>Example</B>
    There exists several constructors for this class one, 
    where you can provide the whole parameter set and one 
    where you can only provide the muscle specific parameters 
    and leave the other values on default value. This example
    shows the version where you only set the muscle specific
    parameters
    @code
    HaeufleTendonDampingCurve fsdeCurve1(1500, 0.12);
    double fsdeValue  = fsdeCurve1.calcValue(serialElasticLengthVelocity, muscleTendonComplexForce);
    @endcode

    <B>References</B>
    \li D.F.B. Haeufle, M. Guenther, A. Bayer, S. Schmitt(2014) Hill-type
        muscle model with serial damping and eccentric force–velocity
        relation. Journal of Biomechanics

    @author Mike Spahr
*/
class OSIMACTUATORS_API HaeufleTendonDampingCurve : public ModelComponent {
    OpenSim_DECLARE_CONCRETE_OBJECT(
            HaeufleTendonDampingCurve, ModelComponent);

public:
    //==============================================================================
    // PROPERTIES
    //==============================================================================
    OpenSim_DECLARE_PROPERTY(max_isometric_force, double,
            "The maximum isometric force of this muscle.");
    OpenSim_DECLARE_PROPERTY(optimal_fiber_length, double,
            "The optimal fiber length of this muscle.");
    OpenSim_DECLARE_PROPERTY(dse_damping_factor, double,
            "The dse damping factor as used in Moerl et al (2012)");
    OpenSim_DECLARE_PROPERTY(rse_damping_factor, double,
            "The rse damping factor as used in Moerl et al (2012)");
    OpenSim_DECLARE_PROPERTY(concentric_contraction_a_rel0, double,
            "The concentric contraction dynamics of CE parameter Arel0");
    OpenSim_DECLARE_PROPERTY(concentric_contraction_b_rel0, double,
            "The concentric contraction dynamics of CE parameter Brel0");

    //==============================================================================
    // PUBLIC METHODS
    //==============================================================================
    /** The default constructor creates a tendon-damping curve using the
    default property values and assigns a default name. */
    HaeufleTendonDampingCurve();

    /** Constructs a tendon-damping curve using the provided muscle specific
     * parameters and assigns a default name. */
    HaeufleTendonDampingCurve(
            double maxIsometricForce, double optimalFiberLength);

    /** Constructs a tendon-damping curve using the provided parameters
    and assigns a default name. */
    HaeufleTendonDampingCurve(double maxIsometricForce,
            double optimalFiberLength, double dseDampingFactor,
            double rseDampingFactor, double concentricContractionArel0,
            double concentricContractionBrel0);

    /** @returns The maximum isometric force of this muscle */
    double getMaxIsometricForce() const;

    /** @returns The optimal fiber length of this muscle */
    double getOptimalFiberLength() const;

    /** @returns The dse damping factor */
    double getDseDampingFactor() const;

    /** @returns The rse damping factor */
    double getRseDampingFactor() const;

    /** @returns The concentric contraction dynamic of CE parameter Arel0 */
    double getConcentricContractionArel0() const;

    /** @returns The concentric contraction dynamic of CE parameter Brel0 */
    double getConcentricContractionBrel0() const;

    /** 
    @param aMaxIsometricForce 
        The maximum isometric force of this muscle        
    */
    void setMaxIsometricForce(double aMaxIsometricForce) const;

    /** 
    @param aOptimalFiberLength
        The optimal fiber length of this muscle
    */
    void setOptimalFiberLength(double aOptimalFiberLength) const;


    /**
    @param aDseDampingFactor
        The Dse damping factor
    @param aRseDampingFactor
        The Rse damping factor
    @param aConcentricContractionArel0
        The concentric contraction dynamic of CE parameter Arel0
    @param aConcentricContractionBrel0
        The concentric contraction dynamic of CE parameter Brel0
    */
    void setTendonDampingParams(double aDseDampingFactor,
            double aRseDampingFactor, double aConcentricContractionArel0,
            double aConcentricContractionBrel0);

    /** Calculates the tendon damping force Fsde which depends on the serial
     * elastic length velocity and the muscleTendonComplexForce */
    double calcValue(double serialElasticLengthVelocity,
            double muscleTendonComplexForce) const;

protected:
    // Component interface.
    void extendFinalizeFromProperties() override;

    //==============================================================================
    // PRIVATE
    //==============================================================================
private:
    void setNull();
    void constructProperties();
};


} // namespace OpenSim

#endif // OPENSIM_HAEUFLE_TENDON_DAMPING_CURVE_H_