/* -------------------------------------------------------------------------- *
 *                      OpenSim:  HaeufleForceVelocityCurve.cpp               *
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
#include "HaeufleForceVelocityCurve.h"
#include <OpenSim/Actuators/HaeufleActiveForceLengthCurve.h>

using namespace OpenSim;
using namespace SimTK;
using namespace std;

//==============================================================================
// CONSTRUCTION
//==============================================================================
// Uses default (compiler-generated) destructor, copy constructor, copy
// assignment operator.
HaeufleForceVelocityCurve::HaeufleForceVelocityCurve() {
    setNull();
    constructProperties();
    setName(getConcreteClassName());
}

HaeufleForceVelocityCurve::HaeufleForceVelocityCurve(
        double concentricContractionARel0, double concentricContractionBRel0,
        double maxForceEccentricExtension, double slopeFactor) 
{
    setNull();
    constructProperties();
    setName(getConcreteClassName());

    set_concentric_contraction_a_rel0(concentricContractionARel0);
    set_concentric_contraction_b_rel0(concentricContractionBRel0);
    set_max_force_eccentric_extension(maxForceEccentricExtension);
    set_slopefactor(slopeFactor);
}

void HaeufleForceVelocityCurve::setNull() { setAuthors("Mike Spahr"); }

void HaeufleForceVelocityCurve::constructProperties() 
{
    constructProperty_concentric_contraction_a_rel0(0.2);
    constructProperty_concentric_contraction_b_rel0(2.0);
    constructProperty_max_force_eccentric_extension(1.5);
    constructProperty_slopefactor(2.0);
}

//==============================================================================
// GET AND SET METHODS
//==============================================================================
double HaeufleForceVelocityCurve::getConcentricContractionARel0() const {
    return get_concentric_contraction_a_rel0();
}

double HaeufleForceVelocityCurve::getConcentricContractionBRel0() const {
    return get_concentric_contraction_b_rel0();
}

double HaeufleForceVelocityCurve::getMaxForceEccentricExtension() const {
    return get_max_force_eccentric_extension();
}

double HaeufleForceVelocityCurve::getSlopeFactor() const {
    return get_slopefactor();
}

void HaeufleForceVelocityCurve::setCurveShape(
    double aConcentricContractionARel0, double aConcentricContractionBRel0, 
    double aMaxForceEccentricExtension, double aSlopeFactor) const 
{
    set_concentric_contraction_a_rel0(aConcentricContractionARel0);
    set_concentric_contraction_b_rel0(aConcentricContractionBRel0);
    set_max_force_eccentric_extension(aMaxForceEccentricExtension);
    set_slopefocator(aSlopeFactor);
}

double HaeufleForceVelocityCurve::calcValueWithoutFmax(double normFiberVelocity,
        double normFiberLength, double activation,
        double correspondingForceLengthValue,
        const HaeufleActiveForceLengthCurve& HaeufleActiveForceLengthCurve)
        const {
    double Fce = 0; // initial Fce to 0
    double Fisom = HaeufleActiveForceLengthCurve.calcValue(normFiberLength);
    // check if this is a eccentric of concentric lengthening
    if (normFiberVelocity <= 0) { // concentric lengthening 
        double Arel = calcArel(
            normFiberLength, activation, HaeufleActiveForceLengthCurve);
        double Brel = calcBrel(activation);
        Fce = (activation * Fisom + Arel) / (1 - normFiberVelocity / Brel) -
              Arel;
    
    } else { // ecentric lengthening
        double Arele = calcArele(activation, correspondingForceLengthValue);
        double Brele = calcBrele(normFiberLength, activation,
                correspondingForceLengthValue, HaeufleActiveForceLengthCurve);
        Fce = (activation * Fisom + Arele) / (1 - normFiberVelocity / Brele) -
              Arele;
    }
    return Fce;
}

double HaeufleForceVelocityCurve::calcArel(double normFiberLength,
        double activation,
        const HaeufleActiveForceLengthCurve& HaeufleActiveForceLengthCurve) const {
    double Arel0 = getConcentricContractionARel0();
    double Qarel = 1 / 4 * (1 + 3 * activation);
    double Larel = 0; // initialize to set total Arel to zero 
    if (normFiberLength < 1) {
        Larel = 1;
    } else {
        Larel = HaeufleActiveForceLengthCurve.calcValue(normFiberLength);
    }
    return Arel0 * Qarel * Larel;
}

double HaeufleForceVelocityCurve::calcBrel(double activation) const 
{
    double Brel0 = getConcentricContractionBRel0();
    double Qbrel = 1 / 7 * (3 + 4 * activation);
    // double Lbrel = 1; // not necessary
    return Brel0 * Qbrel;
}

double HaeufleForceVelocityCurve::calcArele(
    double activation, double correspondingForceLengthValue) const 
{
    double maxEccentricForce = getMaxForceEccentricExtension();
    double Arele =
            -maxEccentricForce * activation * correspondingForceLengthValue;
    return Arele;
}

double HaeufleForceVelocityCurve::calcBrele(double normFiberLength,
        double activation, double correspondingForceLengthValue,
        const HaeufleActiveForceLengthCurve& HaeufleActiveForceLengthCurve) const {
    double Arel = calcArel(normFiberLength, activation,
            HaeufleActiveForceLengthCurve& HaeufleActiveForceLengthCurve);
    double Brel = calcBrel(activation);
    double slopefactor = getSlopeFactor();
    double maxEccentricForce = getMaxForceEccentricExtension();

    double Brele =
            Brel * (1 - maxEccentricForce) /
            (slopefactor * (1 + Arel / (activation * correspondingForceLengthValue)));
    return Brel;
}

