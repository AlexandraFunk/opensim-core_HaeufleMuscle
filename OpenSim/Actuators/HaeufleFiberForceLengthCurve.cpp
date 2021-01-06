/* -------------------------------------------------------------------------- *
 *                    OpenSim:  HaeufleFiberForceLengthCurve.cpp              *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Copyright (c) 2005-2017 Stanford University and the Authors                *
 * Author(s): Mike Spahr	                                                  *
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
#include "HaeufleFiberForceLengthCurve.h"

using namespace OpenSim;
using namespace SimTK;
using namespace std;


//==============================================================================
// CONSTRUCTION
//==============================================================================
// Uses default (compiler-generated) destructor, copy constructor, and copy
// assignment.
HaeufleFiberForceLengthCurve::HaeufleFiberForceLengthCurve() {
    setNull();
    constructProperties();
    setName(getConcreteClassName());
}

HaeufleFiberForceLengthCurve::HaeufleFiberForceLengthCurve(
        double parallelElasticZeroLength, double parallelElasticExponent,
        double parallelElasticForceRelToFmax) {
    setNull();
    constructProperties();
    setName(getConcreteClassName());

    set_parallel_elastic_zero_length(parallelElasticZeroLength);
    set_parallel_elastic_exponent(parallelElasticExponent);
    set_parallel_elastic_force_rel_to_fmax(parallelElasticForceRelToFmax);
}

void HaeufleFiberForceLengthCurve::setNull() { setAuthors("Mike Spahr"); }

void HaeufleFiberForceLengthCurve::constructProperties() {
    constructProperty_parallel_elastic_zero_length(0.95);
    constructProperty_parallel_elastic_exponent(2.5);
    constructProperty_parallel_elastic_force_rel_to_fmax(2.0);
}

void HaeufleFiberForceLengthCurve::extendFinalizeFromProperties() 
{
    Super::extendFinalizeFromProperties();

    std::string errorLocation =
            getName() +
            " HaeufleFiberForceLengthCurve::extendFinalizeFromProperties";

    // Ensure property values are within appropriate ranges.
    OPENSIM_THROW_IF_FRMOBJ(get_parallel_elastic_zero_length() <= 0,
            InvalidPropertyValue,
            getProperty_parallel_elastic_zero_length().getName(),
            "The zero length of the parallel elastic element must be greater "
            "than zero");
    OPENSIM_THROW_IF_FRMOBJ(get_parallel_elastic_exponent() <= 0,
            InvalidPropertyValue,
            getProperty_parallel_elastic_exponent().getName(),
            "The exponent of the parallel elastic element must be greater than "
            "zero");
    OPENSIM_THROW_IF_FRMOBJ(get_parallel_elastic_force_rel_to_fmax() <= 0,
            InvalidPropertyValue,
            getProperty_parallel_elastic_force_rel_to_fmax().getName(),
            "The elastic force relative to the maximum isometric force of the "
            "parallel elastic element must be greater than zero");
}

double HaeufleFiberForceLengthCurve::getParallelElasticZeroLength() const {
    return get_parallel_elastic_zero_length();
}

double HaeufleFiberForceLengthCurve::getParallelElasticExponent() const {
    return get_parallel_elastic_exponent();
}

double HaeufleFiberForceLengthCurve::getParallelElasticForceRelToFmax() const {
    return get_parallel_elastic_force_rel_to_fmax();
}

void HaeufleFiberForceLengthCurve::setParallelElasticParameters(
    double aParallelElasticZeroLength, double aParallelElasticExponent,
    double aParallelElasticForceRelToFmax) {
    set_parallel_elastic_zero_length(aParallelElasticZeroLength);
    set_parallel_elastic_exponent(aParallelElasticExponent);
    set_parallel_elastic_force_rel_to_fmax(aParallelElasticForceRelToFmax);
}

double HaeufleFiberForceLengthCurve::calcValue(
        double FiberLength, double maxIsometricForce, double optimalFiberLength,
        double widthDescendingLimb) const {
    double Lpee0 = getParallelElasticZeroLength() * optimalFiberLength;
    double Fpee = 0; // standard case that FiberLength is smaller than Lpee0
    // check if Fiberlength is larger than or equal to Lpee0
    if (FiberLength >= Lpee0) 
    {
        double Kpee = calculateKPEE(
                maxIsometricForce, optimalFiberLength, widthDescendingLimb);
        Fpee = Kpee * exp((FiberLength - Lpee0), getParallelElasticExponent());  
    } 
    return Fpee;
}

double HaeufleFiberForceLengthCurve::calculateKPEE(double maxIsometricForce,
    double optimalFiberLength, double widthDescendingLimb) {
    double Fpee = getParallelElasticForceRelToFmax();
    double Lpee0factor = getParallelElasticZeroLength();
    double nuepee = getParallelElasticExponent();

    double Kpee = Fpee * maxIsometricForce /
            (exp(optimalFiberLength * (widthDescendingLimb + 1 - Lpee0factor),
                          nuepee));
    return Kpee;
}
