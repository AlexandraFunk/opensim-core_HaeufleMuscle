#ifndef OPENSIM_HAEUFLE_FIBER_FORCE_LENGTH_CURVE_H_
#define OPENSIM_HAEUFLE_FIBER_FORCE_LENGTH_CURVE_H_
/* -------------------------------------------------------------------------- *
 *                     OpenSim:  HaeufleFiberForceLengthCurve.h               *
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
/** This class serves as a FiberForceLengthCurve as described by Haeufle et 
    al in Hill-type muscle model with serial damping and eccentric force–velocity 
    relation. Journal of Biomechanics (D.F.B. Haeufle, M. Guenther, A. Bayer, S. 
    Schmitt (2014)) 
    In contrast to the quintic bezier splines, which is used in the 
    FiberForceLengthCurve, this implementation uses a segmented function.
    And depends on these three free parameters:

    @param ParallelElasticZeroLength
        Zero length of the parallel elastic element as used 
        in Haeufle et al (2014).
    @param ParallelElasticExponent
        Exponent of the parallel elastic elemnt as used in 
        Haeufle et al (2014).
    @param ParallelElasticForceRelToFmax
        Parallel elastic element force relative to Fmax at 
        Lceopt*(1+dWdes) as used in Haeufle et al (2014).

    The required parameters can be set using either the constructor or the
    setParallelElasticParameters function;

    <B>Default Parameter Values</B>
    \verbatim
    ParallelElasticZeroLength ....... 0.95
    ParallelElasticExponent.......... 2.5
    ParallelElasticForceRelToFmax ... 2.0
    \endverbatim

    <B>Example</B>
    @code
    HaeufleFiberForceLengthCurve fpeCurve1;
    fpeCurve1.setParallelElasticParameters(0.95, 2.5, 2.0);
    double fpeVal1 = fpeCurve1.calcValue(0.2);
    @endcode

    <B>References</B>
    \li D.F.B. Haeufle, M. Guenther, A. Bayer, S. Schmitt(2014) Hill-type
        muscle model with serial damping and eccentric force–velocity
        relation. Journal of Biomechanics

    @author Mike Spahr
*/
class OSIMACTUATORS_API HaeufleFiberForceLengthCurve : public ModelComponent {
    OpenSim_DECLARE_CONCRETE_OBJECT(
            HaeufleFiberForceLengthCurve, ModelComponent);

public:
    //==============================================================================
    // PROPERTIES
    //==============================================================================
    // TODO add width of descending fisom branch as well as property??
    OpenSim_DECLARE_PROPERTY(parallel_elastic_zero_length, double,
            "Zero length of the parallel elastic element");
    OpenSim_DECLARE_PROPERTY(parallel_elastic_exponent, double,
            "Exponent of the parallel elastic element");
    OpenSim_DECLARE_OPTIONAL_PROPERTY(parallel_elastic_force_rel_to_fmax,
            double,
            " Parallel elastic element force relative to Fmax at "
            "Lceopt*(1+dWdes)");
    //==============================================================================
    // PUBLIC METHODS
    //==============================================================================
    /** The default constructor creates a fiber-force-length curve using the
    default property values and assigns a default name. */
    HaeufleFiberForceLengthCurve();

    /** Constructs a fiber-force-length curve using the provided parameters and
    assigns a default name. See class documentation for the meaning of these
    parameters, each of which corresponds to a property. */
    HaeufleFiberForceLengthCurve(double parallelElasticZeroLength,
            double parallelElasticExponent,
            double parallelElasticForceRelToFmax);


    /** @returns The parallel elastic elements zero length (L_PEE,0) */
    double getParallelElasticZeroLength() const;
    
    /** @returns The parallel elastic elements exponent (nue_PEE) */
    double getParallelElasticExponent() const;
    
    /** @returns The parallel elastic elemnts force relative to Fmax at 
        Lceopt*(1+dWdes) */
    double getParallelElasticForceRelToFmax() const;

    /**
    @param aParallelElasticZeroLength
        The parallel elastic elements zero lenght
    @param aParallelElasticExponent
        The parallel elastic elements exponent
    @param aParallelElasticForceRelToFmax
        The parallel elastic elemnsts force relative to Fmax at Lceopt*(1+dWdes)
    */
    void setParallelElasticParameters(double aParallelElasticZeroLength,
            double aParallelElasticExponent,
            double aParallelElasticForceRelToFmax);

    /** Evaluates the fiber-force-length curve at a fiber length of
     * 'FiberLength'. */
    double calcValue(double FiberLength, double maxIsometricForce,
            double optimalFiberLength, double widthDescendingLimb) const;

protected:
    // Component interface.
    void extendFinalizeFromProperties() override;

    //==============================================================================
    // PRIVATE
    //==============================================================================
private:
    void setNull();
    void constructProperties();

    /**
    @param maxIsometricForce
        The maximum Isometric force of this muscle
    @param optimalFiberLength
        The optimal fiber length of this muscle
    @param widthDescendingLimb
        The width of the descending limb of the normalized bell curve.
    @returns The constanst K for the parallel elastic element which is 
        neded during the calceValue routine 
    */
    double calculateKPEE(double maxIsometricForce, double optimalFiberLength,
            double widthDescendingLimb);
};

} // namespace OpenSim

#endif // OPENSIM_HAEUFLE_FIBER_FORCE_LENGTH_CURVE_H_