/* -------------------------------------------------------------------------- *
 *                          OpenSim:  WrapViaEllipse.cpp                   *
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


#include "WrapViaEllipse.h"
#include "PathWrap.h"
#include "WrapResult.h"
#include <OpenSim/Common/SimmMacros.h>
#include <OpenSim/Common/Mtx.h>
#include <OpenSim/Common/ModelDisplayHints.h>
#include <OpenSim/Common/ScaleSet.h>

using namespace std;
using namespace OpenSim;
using SimTK::Vec3;

static const char* wrapTypeName = "ViaEllipse";

//=============================================================================
// CONSTRUCTOR(S) AND DESTRUCTOR
//=============================================================================
//_____________________________________________________________________________
/**
 * Default constructor.
 */
WrapViaEllipse::WrapViaEllipse()
{ 
	constructProperties();
}

//_____________________________________________________________________________
/**
 * Destructor.
 */
WrapViaEllipse::~WrapViaEllipse()
{

}

//_____________________________________________________________________________
/**
 * Connect properties to local pointers.
 */
void WrapViaEllipse::constructProperties() 
{
    constructProperty_semi_axis_length_H(0.05); // TODO set realistic default value
    constructProperty_semi_axis_length_G(0.05); // TODO set realistic default value
}

void WrapViaEllipse::extendFinalizeFromProperties()
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

void WrapViaEllipse::extendAddToSystem(SimTK::MultibodySystem& system) const 
{
    Super::extendAddToSystem(system);

    this->_viaEllipseInfoCV = addCacheVariable("viaEllipsePlottingInfos",
            ViaEllipsePlottingInfos(), SimTK::Stage::Velocity);
}

const double WrapViaEllipse::getSemiAxisLengthH() const 
{
    return get_semi_axis_length_H();
}

void WrapViaEllipse::setSemiAxisLengthH(double aSemiAxisLengthH) 
{
    set_semi_axis_length_H(aSemiAxisLengthH);
}

const double WrapViaEllipse::getSemiAxisLengthG() const 
{
    return get_semi_axis_length_G();
}

void WrapViaEllipse::setSemiAxisLengthG(double aSemiAxisLengthG) 
{
    set_semi_axis_length_G(aSemiAxisLengthG);
}

double WrapViaEllipse::getAngleOnEllipse(const SimTK::State& s) const {
    return getViaEllipsePlottingInfos(s).phi;
}

double WrapViaEllipse::getDeflectionPointX(const SimTK::State& s) const {
    return getViaEllipsePlottingInfos(s).deflectionPoint[0];
}

double WrapViaEllipse::getDeflectionPointY(const SimTK::State& s) const {
    return getViaEllipsePlottingInfos(s).deflectionPoint[1];
}

double WrapViaEllipse::getDeflectionPointZ(const SimTK::State& s) const {
    return getViaEllipsePlottingInfos(s).deflectionPoint[2];
}

const WrapViaEllipse::ViaEllipsePlottingInfos&
WrapViaEllipse::getViaEllipsePlottingInfos(const SimTK::State& s) const {
    if (isCacheVariableValid(s, _viaEllipseInfoCV)) {
        return getCacheVariableValue(s, _viaEllipseInfoCV);
    }

    /* not valid variable */
    ViaEllipsePlottingInfos& umdi = updCacheVariableValue(s, _viaEllipseInfoCV);
    //calcMuscleDynamicsInfo(s, umdi);
    //markCacheVariableValid(s, _viaEllipseInfoCV);
    return umdi;
}


void WrapViaEllipse::extendScale(
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
const char* WrapViaEllipse::getWrapTypeName() const 
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
std::string WrapViaEllipse::getDimensionsString() const 
{
    stringstream dimensions;
    dimensions << "semi axis length H " << get_semi_axis_length_H()
               << "\nsemi axis length H: " << get_semi_axis_length_G();
    return dimensions.str();
}

/// Implement generateDecorations to draw geometry in visualizer
void WrapViaEllipse::generateDecorations(bool fixed, const ModelDisplayHints& hints,
    const SimTK::State& state, SimTK::Array_<SimTK::DecorativeGeometry>& appendToThis) const 
{
    Super::generateDecorations(fixed, hints, state, appendToThis);
    if (!fixed) 
    {
        return;
    }

    if (hints.get_show_wrap_geometry()) {
        const Appearance& defaultAppearance = get_Appearance();
        if (!defaultAppearance.get_visible()) return;
        const Vec3 color = defaultAppearance.get_color();

        /** 
        * since DecorativeEllipsoid function generates a ellipsoid around the 
        * center of the ellipse frame, which is basically defied as the rattachment
        * point. Therefore, we need to shift the ellipse frame to the real center 
        * of the ellipse, which lays the amount of the semi axis length H in 
        * positive y-direction.
        **/ 
        SimTK::Transform translation_to_ell_center(SimTK::Vec3(0., get_semi_axis_length_H(), 0.));
        const SimTK::Transform X_BP = calcWrapGeometryTransformInBaseFrame() *
                                      translation_to_ell_center;
        /**
        * create decorative Ellipsoid since there is no decorative ellipse
        * available in SimTK. Here we just assume that the x dimension of the
        * ellipsoid is very small, making it look similar to an 2D object
        * Since this is only for visualization, the visual representation will
        * not change the behaviour of the algorithm.
        **/
        const Vec3 ellipsoidShape =
                Vec3(0.001, getSemiAxisLengthH(), getSemiAxisLengthG());

        appendToThis.push_back(
            SimTK::DecorativeEllipsoid(ellipsoidShape)
            .setTransform(X_BP).setResolution(2.0)
            .setColor(color).setOpacity(defaultAppearance.get_opacity())
            .setScale(1).setRepresentation(defaultAppearance.get_representation())
            .setBodyId(getFrame().getMobilizedBodyIndex()));
    }
}

int WrapViaEllipse::wrapLine(const SimTK::State& s, SimTK::Vec3& aPoint1,
        SimTK::Vec3& aPoint2,
        const PathWrap& aPathWrap, WrapResult& aWrapResult, bool& aFlag) const 
{
    double tmp, theta;
    SimTK::Vec3 e_x, e_y, e_z, M, M_defl(0., 0., 0.), G_defl(0., 0., 0.),
            H_defl(0., 0., 0.);
    SimTK::Vec3 placeholder;
    double phi = 0;

    SimTK::Vec3 H(0., get_semi_axis_length_H(), 0.);
    SimTK::Vec3 G(0., 0., -get_semi_axis_length_G());


    // start with calculation point 1 and point 2 into ellipse deflection frame for easier calculation
    // aPoint1 and a Point2 -> defined in the ellipse frame
    // semi lengths g and h -> defined as scalars in the ellipse frame
    M = H - aPoint1;

    Mtx::Normalize(3, aPoint2 - aPoint1, e_x);
    Mtx::CrossProduct(aPoint2 - aPoint1, M, e_z);

    // if M is a point on line aPoint2aPoint1, ellipse is neglected
    if (e_z == 0) { 
        aFlag = false;
        return 0; 
    }

    Mtx::Normalize(3, e_z, e_z);
    Mtx::CrossProduct(e_z, e_x, e_y);

    // calculate new center of deflection frame and express G and H in deflection frame
    M_defl[0] = Mtx::DotProduct(3, M, e_x);
    M_defl[1] = Mtx::DotProduct(3, M, e_y);
    G_defl[0] = Mtx::DotProduct(3, G, e_x);
    G_defl[1] = Mtx::DotProduct(3, G, e_y);
    G_defl[2] = Mtx::DotProduct(3, G, e_z);
    H_defl[0] = Mtx::DotProduct(3, H, e_x);
    H_defl[1] = Mtx::DotProduct(3, H, e_y);
    H_defl[2] = Mtx::DotProduct(3, H, e_z);

    // check if ellipse can be neglected or not
    if (abs(H_defl[2]) < SimTK::SignificantReal) {
        tmp = M_defl[1] / H_defl[1];
    } else {
        if (abs(G_defl[2]) < SimTK::SignificantReal) {
            tmp = -M_defl[1] / G_defl[1];
        } else {
            theta = atan(H_defl[2] / G_defl[2]);
            tmp = -M_defl[1] / sin(theta) /
                  (G_defl[1] - H_defl[1] * G_defl[2] / H_defl[2]);
        }
    }

    if (abs(tmp) < 1.0) { 
        aFlag = false;
        ViaEllipsePlottingInfos& myInfos =
                updCacheVariableValue(s, _viaEllipseInfoCV);
        // populate struct into Plotting Infos
        myInfos.phi = phi;
        myInfos.deflectionPoint = (0,0,0);
        markCacheVariableValid(s, _viaEllipseInfoCV);
        return 0;
    }
    // else not necessary since from here ellipse is not neglected

    double G_norm = Mtx::Normalize(3, G_defl, placeholder);
    double H_norm = Mtx::Normalize(3, H_defl, placeholder);
    double error = 2.0e-3;
    double s01 = Mtx::Normalize(3, aPoint2 - aPoint1, placeholder);
    SimTK::Vec3 P_defl(0., 0., 0.);

    // first check if ellipse is basically a circle -> calculation is
    // much easier
    if (abs(G_norm - H_norm) < error) {
        if (H_defl[1] == 0) { 
            double phi1 = SimTK::Pi / 2.0;
            double phi2 = 3 * SimTK::Pi / 2.0;
            SimTK::Vec3 P1 =
                    M_defl + G_defl * sin(phi1) - H_defl * cos(phi1);
            SimTK::Vec3 P2 =
                    M_defl + G_defl * sin(phi2) - H_defl * cos(phi2);
            double len1 = pathLengthTroughEllipse(s01, P1);
            double len2 = pathLengthTroughEllipse(s01, P2);
            phi = (len1 > len2) ? phi2 : phi1;
        } else {
            double phi1 = atan(-G_defl[1] / H_defl[1]);
            double phi2 = atan(-G_defl[1] / H_defl[1]) + SimTK::Pi;
            SimTK::Vec3 P1 =
                    M_defl + G_defl * sin(phi1) - H_defl * cos(phi1);
            SimTK::Vec3 P2 =
                    M_defl + G_defl * sin(phi2) - H_defl * cos(phi2);
            double len1 = pathLengthTroughEllipse(s01, P1);
            double len2 = pathLengthTroughEllipse(s01, P2);
            phi = (len1 > len2) ? phi2 : phi1;
        }
    } else {
        // From here numerically calculation of shortest path
        // find starting point on ellipse (phi=0,pi,pi/2 or 3*pi/2)
        P_defl = M_defl + G_defl * sin(0.0) - H_defl * cos(0.0);
        double min_length = pathLengthTroughEllipse(s01, P_defl);
        double len_h1 = 0.0, len_h2 = 0.0, Alpha = SimTK::Pi / 2;

        for (int i = 1; i < 4; i++) {
            P_defl = M_defl + G_defl * sin(SimTK::Pi / 2 * i) -
                H_defl * cos(SimTK::Pi / 2 * i);
            len_h1 = pathLengthTroughEllipse(s01, P_defl);
            if (len_h1 < min_length) {
                min_length = len_h1;
                phi = SimTK::Pi / 2 * i;
            }
        }
        while (Alpha > 0.005) { // error of phi > 0.005 
            P_defl = M_defl + G_defl * sin(phi - Alpha) - H_defl * cos(phi - Alpha);
            len_h1 = pathLengthTroughEllipse(s01, P_defl);
            P_defl = M_defl + G_defl * sin(phi + Alpha) - H_defl * cos(phi + Alpha);
            len_h2 = pathLengthTroughEllipse(s01, P_defl);

            if (len_h1 < min_length) { 
                phi -= Alpha;
                min_length = len_h1;
            }
            if (len_h2 < min_length) {
                phi += Alpha;
                min_length = len_h2;
            }
            Alpha /= 2;
        }            
    }

    // shortest point on ellipse surface
    SimTK::Vec3 pointOnEllipse = H + G * sin(phi) - H * cos(phi); 
   
    ViaEllipsePlottingInfos& myInfos =
            updCacheVariableValue(s, _viaEllipseInfoCV);
    // populate struct into Plotting Infos
    myInfos.phi = phi;
    myInfos.deflectionPoint = pointOnEllipse;
    markCacheVariableValid(s, _viaEllipseInfoCV);

    // recalculate P_defl into ellipse frame and add it to wrap_pts array
    aWrapResult.wrap_pts.insert(0, pointOnEllipse);

    // fill r1 and r2 of WrapResult with NaN value to check if they are used
    aWrapResult.r1 = pointOnEllipse;
    aWrapResult.r2 = SimTK::Vec3(SimTK::NaN, SimTK::NaN, SimTK::NaN);

    aFlag = true;
    return 1;
}

// calculates the path length with one deflecting ellipse
double WrapViaEllipse::pathLengthTroughEllipse(
    double s01, SimTK::Vec3 P) const {
    SimTK::Vec3 placeholder;
    double p_norm = Mtx::Normalize(3, P, placeholder);
    return p_norm + sqrt(pow(P[0] - s01, 2) + pow(P[1], 2) + pow(P[2], 2));
}
