/* -------------------------------------------------------------------------- *
 *                          OpenSim:  WrapHammerEllipse.cpp                   *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Copyright (c) 2005-2017 Stanford University and the Authors                *
 * Author(s): Maria Hammer and Mike Spahr                                     *
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


#include "WrapHammerEllipse.h"

using namespace std;
using namespace OpenSim;
using SimTK::Vec3;

static const char* wrapTypeName = "HammerEllipse";

//=============================================================================
// CONSTRUCTOR(S) AND DESTRUCTOR
//=============================================================================
//_____________________________________________________________________________
/**
 * Default constructor.
 */
WrapHammerEllipse::WrapHammerEllipse()
{ 
	constructProperties();
}

//_____________________________________________________________________________
/**
 * Destructor.
 */
WrapHammerEllipse::~WrapHammerEllipse()
{

}

//_____________________________________________________________________________
/**
 * Connect properties to local pointers.
 */
void WrapHammerEllipse::constructProperties() 
{
    constructProperty_semi_axis_length_H(0.05); // TODO set realistic default value
    constructProperty_semi_axis_length_G(0.05); // TODO set realistic default value
}

void WrapHammerEllipse::extendFinalizeFromProperties()
{
    // Base class
    Super::extendFinalizeFromProperties();

    OPENSIM_THROW_IF_FRMOBJ(get_semi_axis_length_H() < 0, InvalidPropertyValue,
            getProperty_semi_axis_length_H().getName(),
            "Semi axis length H cannot be less than zero");
    OPENSIM_THROW_IF_FRMOBJ(get_semi_axis_length_G() < 0, InvalidPropertyValue,
            getProperty_semi_axis_length_H().getName(),
            "Semi axis length G cannot be less than zero");
}

const double WrapHammerEllipse::getSemiAxisLengthH() const 
{
    return get_semi_axis_length_H();
}

void WrapHammerEllipse::setSemiAxisLengthH(double aSemiAxisLengthH) 
{
    set_semi_axis_length_H(aSemiAxisLengthH);
}

const double WrapHammerEllipse::getSemiAxisLengthG() const 
{
    return get_semi_axis_length_G();
}

void WrapHammerEllipse::setSemiAxisLengthG(double aSemiAxisLengthG) 
{
    set_semi_axis_length_G(aSemiAxisLengthG);
}


void WrapHammerEllipse::extendScale(
    const SimTK::State& s, const ScaleSet& scaleSet) 
{
    Super::extendScale(s, scaleSet);

    // Get scale factors (if an entry for the Frame's base Body exists).
    const Vec3& scaleFactors = getScaleFactors(scaleSet, getFrame());
    if (scaleFactors == ModelComponent::InvalidScaleFactors) 
    { 
        return; 
    }
    
    // _pose.x() holds the ellipse's X-axis expressed in the body's reference
    // frame. The elementwise product of this vector and the scaleFactors vector
    // gives the amount that the ellipsoid must be scaled in the X dimension.
    // Similar for the Y and Z dimensions.
    Vec3 localScaleVector[3];

    localScaleVector[0] = _pose.x().elementwiseMultiply(scaleFactors);
    localScaleVector[1] = _pose.y().elementwiseMultiply(scaleFactors);
    localScaleVector[2] = _pose.z().elementwiseMultiply(scaleFactors);

    SimTK::Vec3 previousDimensions(0, getSemiAxisLengthG(), getSemiAxisLengthH());
    for (int i = 0; i < 3; ++i) 
    {
        previousDimensions[i] *= localScaleVector[i].norm();
    }

    setSemiAxisLengthG(previousDimensions[1]);
    setSemiAxisLengthH(previousDimensions[2]);
}

//_____________________________________________________________________________
/**
 * Get the name of the type of wrap object ("cylinderObst" in this case)
 *
 * @return A string representing the type of wrap object
 */
const char* WrapHammerEllipse::getWrapTypeName() const 
{ 
    return wrapTypeName; 
}

//_____________________________________________________________________________
/**
 * Get a string holding the dimensions definition that SIMM would
 * use to describe this object. This is a rather ugly convenience
 * function for outputting SIMM joint files.
 *
 * @return A string containing the dimensions of the wrap object
 */
std::string WrapHammerEllipse::getDimensionsString() const 
{
    stringstream dimensions;
    dimensions << "semi axis length H " << get_semi_axis_length_H()
               << "\nsemi axis length H: " << get_semi_axis_length_G();
    return dimensions.str();
}

/// Implement generateDecorations to draw geometry in visualizer
void WrapHammerEllipse::generateDecorations(bool fixed,
        const ModelDisplayHints& hints,
    const SimTK::State& state,
    SimTK::Array_<SimTK::DecorativeGeometry>& appendToThis) const 
{
    Super::generateDecorations(fixed, hints, state, appendToThis);
    if (!fixed) return;

    if (hints.get_show_wrap_geometry()) {
        const Appearance& defaultAppearance = get_Appearance();
        if (!defaultAppearance.get_visible()) return;
        const Vec3 color = defaultAppearance.get_color();

        const auto X_BP = calcWrapGeometryTransformInBaseFrame();
        /**
        * create decorative Ellipsoid since there is no decorative ellipse
        * available in SimTK. Here we just assume that the x dimension of the
        * ellipsoid is very small, making it look similar to an 2D object
        * Since this is only for visualization, the visual representation will
        * not change the behaviour of the algorithm.
        **/
        const Vec3 ellipsoidShape =
                Vec3(0.01, getSemiAxisLengthG(), getSemiAxisLengthH());
        appendToThis.push_back(SimTK::DecorativeEllipsoid(ellipsoidShape))
            .setTransform(X_BP).setResolution(2.0)
            .setColor(color).setOpacity(defaultAppearance.get_opacity())
            .setScale(1).setRepresentation(defaultAppearance.get_representation())
            .setBodyId(getFrame().getMobilizedBodyIndex()));
    }
}

int WrapHammerEllipse::wrapLine(const SimTK::State& s, SimTK::Vec3& aPoint1,
        SimTK::Vec3& aPoint2,
        const PathWrap& aPathWrap, WrapResult& aWrapResult, bool& aFlag) const {
    return 0;
}

void WrapHammerEllipse::connectToModelAndBody(
        Model& aModel, PhysicalFrame& aBody) 
{

}
