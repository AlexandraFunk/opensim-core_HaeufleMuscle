#ifndef OPENSIM_HAEUFLE_FORCE_VELOCITY_CURVE_H_
#define OPENSIM_HAEUFLE_FORCE_VELOCITY_CURVE_H_

/* -------------------------------------------------------------------------- *
 *                       OpenSim:  HaeufleForceVelocityCurve.h                *
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
/** This class serves as a ForceVelocityCurve as described by Haeufle et 
    al in Hill-type muscle model with serial damping and eccentric force–velocity 
    relation. Journal of Biomechanics (D.F.B. Haeufle, M. Guenther, A. Bayer, S. 
    Schmitt (2014)) 
    In contrast to the quintic bezier splines, which is used in the 
    ForceVelocityCurve, this implementation uses two hyperbolas to model the fiber
    contraction velocities for concentric and eccentric contractions.
    The provided fiber length must be normalized with the optimalFiberLength.
    The Haeufle force-velocity curve is constructed from 4 properties:
    
    @param concentricContractionARel0
    Derived from the classical Hill constant a

    @param concentricContractionBRel0
    Derived from the classical Hill constant b

    @param maxForceEccentricExtension
    Maximal Force during eccentric lengthening

    @param slopeFactor
    The ratio between the eccentric and concentric derivatives dF/dVce

    <B>Default Parameter Values</B>
    \verbatim
    concentricContractionARel0 .... 0.2
    concentricContractionBRel0 .... 2.0
    maxForceEccentricExtension .... 1.5
    slopeFactor ................... 2.0
    \endverbatim

    <B>Note</B>
    This curve doesn't scale the muscle with the maximum isometric force!

    <B>Example</B>
    @code
        HaeufleForceVelocityCurve fvCurve(0.0, 0.25, 5.0, 0.0);
        double falVal  = fvCurve.calcValue(1.0);
        double dfalVal = fvCurve.calcDerivative(1.0, 1);
    @endcode

    <B>References</B>
    \li D.F.B. Haeufle, M. Guenther, A. Bayer, S. Schmitt(2014) Hill-type
        muscle model with serial damping and eccentric force–velocity
        relation. Journal of Biomechanics

    @author Mike Spahr
*/
class OSIMACTUATORS_API HaeufleForceVelocityCurve : public ModelComponent {
    OpenSim_DECLARE_CONCRETE_OBJECT(HaeufleForceVelocityCurve, ModelComponent);

public:
    //==============================================================================
    // PROPERTIES
    //==============================================================================
    OpenSim_DECLARE_PROPERTY(concentric_contraction_a_rel0, double,
            "Derived from the classical Hill constant a");
    OpenSim_DECLARE_PROPERTY(concentric_contraction_b_rel0, double,
            "Derived from the classical Hill constant b");
    OpenSim_DECLARE_PROPERTY(max_force_eccentric_extension, double, 
            "Maximal Force during eccentric lengthening");
    OpenSim_DECLARE_PROPERTY(slopefactor, double, 
        "The ratio between the eccentric and concentric derivatives dF/dVce");

    //==============================================================================
    // PUBLIC METHODS
    //==============================================================================
    /** The default constructor creates a force-velocity curve using the default
    property values and assigns a default name. */
    HaeufleForceVelocityCurve();

    /** Constructs a force-velocity curve using the provided parameters and
    assigns a default name. */
    HaeufleForceVelocityCurve(double concentricContractionARel0, 
        double concentricContractionBRel0, double maxForceEccentricExtension, 
        double slopeFactor);

    /** @returns The concentric contraction parameter which is called Arel in Haeufle 
    et al */
    double getConcentricContractionARel0() const;

    /** @returns The concentric contraction parameter which is called Brel in
    Haeufle et al */
    double getConcentricContractionBRel0() const;

    /** @returns The maximum force at eccentric lengthening which is called Fe in
    Haeufle et al */
    double getMaxForceEccentricExtension() const;
    
    /** @returns The slopefactor of the curve which is called Se in
    Haeufle et al */    
    double getSlopeFactor() const;

    void setCurveShape(double aConcentricContractionARel0,
            double aConcentricContractionBRel0,
            double aMaxForceEccentricExtension, double aSlopeFactor) const;

    /** Evaluates the force-velocity curve at a normalized fiber velocity of
    'normFiberVelocity'. */
    double calcValueWithoutFmax(double normFiberVelocity,
            double normFiberLength, double activation, 
            double correspondingForceLengthValue,
            const HaeufleActiveForceLengthCurve& HaeufleActiveForceLengthCurve)
            const;

    /** Calculate the fiber length and activation dependent normalized Hill
    parameter Arel for the conccentric case */
    double calcArel(double normFiberLength, double activation,
            const HaeufleActiveForceLengthCurve& HaeufleActiveForceLengthCurve)
            const;

    /** Calculate the activation dependent normalized Hill
    parameter Brel for the conccentric case */
    double calcBrel(double activation) const;

    /** Calculate the activation and correspondingForceLengthValue dependent
    normalized Hill parameter Arele for the eccentric case */
    double calcArele(double activation, double correspondingForceLengthValue) const;

    /** Calculate the activation and correspondingForceLengthValue dependent 
    normalized Hill parameter Brele for the eccentric case */
    double calcBrele(double normFiberLength, double activation,
            double correspondingForceLengthValue,
            const HaeufleActiveForceLengthCurve& HaeufleActiveForceLengthCurve)
            const;
   
    //==============================================================================
    // PRIVATE
    //==============================================================================
private:
    void setNull();
    void constructProperties();
};

} // namespace OpenSim

#    endif // OPENSIM_HAEUFLE_FORCE_VELOCITY_CURVE_H_
