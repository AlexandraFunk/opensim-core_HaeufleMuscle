#ifndef OPENSIM_HAEUFLE_TENDON_FORCE_LENGTH_CURVE_H_
#define OPENSIM_HAEUFLE_TENDON_FORCE_LENGTH_CURVE_H_
/* -------------------------------------------------------------------------- *
 *                     OpenSim:  HaeufleTendonForceLengthCurve.h                     *
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
#include <OpenSim/Actuators/osimActuatorsDLL.h>
#include <OpenSim/Simulation/Model/ModelComponent.h>

#ifdef SWIG
    #ifdef OSIMACTUATORS_API
        #undef OSIMACTUATORS_API
        #define OSIMACTUATORS_API
    #endif
#endif

namespace OpenSim {
/** This class serves as a TendonForceLengthCurve as described by Guenther et al
    in 2007. 
    In contrast to the quintic bezier splines, which is used in the 
    TendonForceLengthCurve, this implementation uses a segmentwise defined
    function with the following parameters.

    @param serialElasticRestLength
    @param relativeStretchAtNonlinearLinearTransition
    @param relativeStretchAtLinearPart
    @param forceAtNonlinearLinearTransition
   

    <B>Default Parameter Values</B>
    \verbatim
    serialElasticRestLength ...................... 0.15
    relativeStretchAtNonlinearLinearTransition ... 0.0425
    relativeStretchAtLinearPart .................. 0.0170
    forceAtNonlinearLinearTransition ............. 500
    \endverbatim

    <B>Example</B>
    @code
        HaeufleTendonForceLengthCurve fseCurve(0.15, 0.0425, 0.0170, 500);
        double fseVal   = fseCurve.calcValue(0.02);
    @endcode

    <B>References</B>
    \li D.F.B. Haeufle, M. Guenther, A. Bayer, S. Schmitt(2014) Hill-type
        muscle model with serial damping and eccentric force–velocity
        relation. Journal of Biomechanics

    @author Mike Spahr
*/
class OSIMACTUATORS_API HaeufleTendonForceLengthCurve : public ModelComponent {
    OpenSim_DECLARE_CONCRETE_OBJECT(
            HaeufleTendonForceLengthCurve, ModelComponent);

public:
    //==============================================================================
    // PROPERTIES
    //==============================================================================
    OpenSim_DECLARE_PROPERTY(serial_elastic_rest_length, double,"");
    OpenSim_DECLARE_OPTIONAL_PROPERTY(
            relative_stretch_at_nonlinear_linear_transition, double,"");
    OpenSim_DECLARE_OPTIONAL_PROPERTY(relative_stretch_at_linear_part, double,
            "");
    OpenSim_DECLARE_OPTIONAL_PROPERTY(force_at_nonlinear_linear_transition, double,
            "");

    //==============================================================================
    // PUBLIC METHODS
    //==============================================================================
    /** The default constructor creates a tendon-force-length curve using the
    default property values and assigns a default name. */
    HaeufleTendonForceLengthCurve();

    /** Constructs a tendon-force-length curve using the muscle dependen parameters and
    assigns a default name. */
    HaeufleTendonForceLengthCurve(double serialElasticRestLength,
            double forceAtNonlinearLinearTransition);

    /** Constructs a tendon-force-length curve using the provided parameters and
    assigns a default name. */
    HaeufleTendonForceLengthCurve(double serialElasticRestLength,
            double relativeStretchAtNonlinearLinearTransition,
            double relativeStretchAtLinearPart,
            double forceAtNonlinearLinearTransition);

    /** @returns The serial elastic rest length Lsee,0 */
    double getSerialElasticRestLength() const;
    
    /** @returns The relative stretch at nonlinear/linear transition in Fsee */
    double getRelativeStretchAtNonlinearLinearTransition() const;
    
    /** @returns The relative stretch in the linear part for force increase
     * delta Fsee,0 */
    double getRelativeStretchAtLinearPart() const;
    
    /** @returns The force at nonlinear/linear transition in Fsee */
    double getForceAtNonlinearLinearTransition() const;

    /** 
    @param aSerialElasticRestLength
        The resting length of the serial elastic element
    */
    void setSerialElasticRestLength(double aSerialElasticRestLength ) const;
   
    /**
    @param aForceAtNonlinearLinearTransition
        The force at non-linear/linear transition which is muscle specific
    */
    void setForceAtNonlinearLinearTransition(double aForceAtNonlinearLinearTransition ) const;
    
    /**
    @param aForceAtNonlinearLinearTransition
        The force at non-linear/linear transition which is muscle specific
    */
    void setRelativeStretchAtNonlinearLinearTransition(double aRelativeStretchAtNonlinearLinearTransition ) const;
    
    /**
    @param aForceAtNonlinearLinearTransition
        The force at non-linear/linear transition which is muscle specific
    */
    void setRelativeStretchAtLinearPart(double aRelativeStretchAtLinearPart ) const;


    /** Evaluates the tendon-force-length curve at a tendon length of
    'aSerialElasticLength'. */
    double calcValue(double aSerialElasticLength) const;

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

#endif // OPENSIM_HAEUFLE_TENDON_FORCE_LENGTH_CURVE_H_
