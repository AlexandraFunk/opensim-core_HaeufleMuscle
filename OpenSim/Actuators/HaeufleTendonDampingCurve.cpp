/* -------------------------------------------------------------------------- *
 *                      OpenSim:  HaeufleTendonDampingCurve.cpp               *
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
#include "HaeufleTendonDampingCurve.h"

using namespace OpenSim;
using namespace SimTK;
using namespace std;

//==============================================================================
// CONSTRUCTION
//==============================================================================

HaeufleTendonDampingCurve::HaeufleTendonDampingCurve() 
{ 
    setNull();
    constructProperties();
    setName(getConcreteClassName());
}

HaeufleTendonDampingCurve::HaeufleTendonDampingCurve(
        double maxIsometricForce, double optimalFiberLength) 
{
    setNull();
    constructProperties();
    setName(getConcreteClassName());

    set_max_isometric_force(maxIsometricForce);
    set_optimal_fiber_length(optimalFiberLength);
}

HaeufleTendonDampingCurve::HaeufleTendonDampingCurve(double maxIsometricForce,
        double optimalFiberLength,
        double dseDampingFactor, double rseDampingFactor, double concentricContractionArel0,
        double concentricContractionBrel0) 
{
    setNull();
    constructProperties();
    setName(getConcreteClassName());

    set_max_isometric_force(maxIsometricForce);
    set_optimal_fiber_length(optimalFiberLength);
    set_dse_damping_factor(dseDampingFactor);
    set_rse_damping_factor(rseDampingFactor);
    set_concentric_contraction_a_rel0(concentricContractionArel0);
    set_concentric_contraction_b_rel0(concentricContractionBrel0);
}

void HaeufleTendonDampingCurve::setNull() { setAuthors("Mike Spahr"); }

void HaeufleTendonDampingCurve::constructProperties() {
    constructProperty_max_isometric_force(1500);
    constructProperty_optimal_fiber_length(0.12);
    constructProperty_dse_damping_factor(0.3);
    constructProperty_rse_damping_factor(0.01);
    constructProperty_concentric_contraction_a_rel0(0.2);
    constructProperty_concentric_contraction_b_rel0(2.0);
}

void HaeufleTendonDampingCurve::extendFinalizeFromProperties() {
    Super::extendFinalizeFromProperties();

    std::string errorLocation =
            getName() +
            " HaeufleTendonDampingCurve::extendFinalizeFromProperties";

    // Ensure property values are within appropriate ranges.
    OPENSIM_THROW_IF_FRMOBJ(get_max_isometric_force() <= 0,
            InvalidPropertyValue, getProperty_max_isometric_force().getName(),
            "The maximum isometric force must be greater than zero");
    OPENSIM_THROW_IF_FRMOBJ(get_optimal_fiber_length() <= 0,
            InvalidPropertyValue, getProperty_optimal_fiber_length.getName(),
            "The optimal fiber length must be greater than zero");
    OPENSIM_THROW_IF_FRMOBJ(get_dse_damping_factor() <= 0,
            InvalidPropertyValue,
            getProperty_dse_damping_factor().getName(),
            "The dse damping factor must be greater than zero")
    OPENSIM_THROW_IF_FRMOBJ(get_rse_damping_factor() <= 0, InvalidPropertyValue,
            getProperty_rse_damping_factor().getName(),
            "The rse damping factor must be greater than zero");
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
}

double HaeufleTendonDampingCurve::getMaxIsometricForce() const {
    return get_max_isometric_force();
}

double HaeufleTendonDampingCurve::getOptimalFiberLength() const {
    return get_optimal_fiber_length();
}

double HaeufleTendonDampingCurve::getDseDampingFactor() const {
    return get_dse_damping_factor();
}

double HaeufleTendonDampingCurve::getRseDampingFactor() const {
    return get_rse_damping_factor();
}

double HaeufleTendonDampingCurve::getConcentricContractionArel0() const {
    return get_concentric_contraction_a_rel0();
}

double HaeufleTendonDampingCurve::getConcentricContractionBrel0() const {
    return get_concentric_contraction_b_rel0();
}

void HaeufleTendonDampingCurve::setMaxIsometricForce(
    double aMaxIsometricForce) const {
    set_max_isometric_force(aMaxIsometricForce);
}

void HaeufleTendonDampingCurve::setOptimalFiberLength(
    double aOptimalFiberLength) const {
    set_optimal_fiber_length(aOptimalFiberLength);
}

void HaeufleTendonDampingCurve::setTendonDampingParams(double aDseDampingFactor,
    double aRseDampingFactor,
    double aConcentricContractionArel0, double aConcentricContractionBrel0) {
    set_dse_damping_factor(aDseDampingFactor);
    set_rse_damping_factor(aRseDampingFactor);
    set_concentric_contraction_a_rel0(aConcentricContractionArel0);
    set_concentric_contraction_b_rel0(aConcentricContractionBrel0);
}

double HaeufleTendonDampingCurve::calcValue(double serialElasticLengthVelocity,
    double muscleTendonComplexForce) const {

    double Dse = getDseDampingFactor();
    double Rse = getRseDampingFactor();
    double Fmax = getMaxIsometricForce();
    double lceopt = getOptimalFiberLength();
    double Arel0 = getConcentricContractionArel0();
    double Brel0 = getConcentricContractionBrel0();

    double Fsde = Dse * Fmax * Arel0 / (lceopt * Brel0) *
                  ((1 - Rse) * muscleTendonComplexForce / Fmax + Rse) *
                  serialElasticLengthVelocity;
    return Fsde;
}
