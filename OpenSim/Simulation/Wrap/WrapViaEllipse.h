#ifndef OPENSIM_WRAP_HAMMER_ELLIPSE_H_
#define OPENSIM_WRAP_HAMMER_ELLIPSE_H_
/* -------------------------------------------------------------------------- *
 *                           OpenSim:  WrapViaEllipse.h                    *
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

// INCLUDE
# include "WrapObject.h"

namespace OpenSim {

class PathWrap;
class WrapResult;
//=============================================================================
//=============================================================================
/**
 * A class implementing an empty ellipse obstacle for muscle wrapping, based on
 * algorithm presented in Hammer et al (2019).
 *
 * @author Maria Hammer, Mike Spahr
 */
class OSIMSIMULATION_API WrapViaEllipse : public WrapObject {
    OpenSim_DECLARE_CONCRETE_OBJECT(WrapViaEllipse, WrapObject);

public:
//==============================================================================
// PROPERTIES
//==============================================================================

    /**
    * Schon in WrapObject enthalten??
    *
    *
    OpenSim_DECLARE_PROPERTY(xyz_ellipse_frame_rotation, SimTK::Vec3,
            "Body-fixed Euler angle sequence describing the orientation of the "
            "ellipse reference frame relative to the body reference system.");

    OpenSim_DECLARE_PROPERTY(ellipse_attachment_point, SimTK::Vec3,
        "Location of the ellipse attachment point.");
    **/

    OpenSim_DECLARE_PROPERTY(semi_axis_length_H, double,
            "Ellipse semi-axis length in z-Direction of the ellipse "
            "reference frame.");
    OpenSim_DECLARE_PROPERTY(semi_axis_length_G, double,
            "Ellipse semi-axis length in y-Direction of the ellipse reference "
            "frame.");

public:
//=============================================================================
// METHODS
//=============================================================================

    // constructor
    WrapViaEllipse();
    // destructor
    virtual ~WrapViaEllipse();

    /** Scale the torus's dimensions. The base class (WrapObject) scales the
    origin of the torus in the body's reference frame. */
    void extendScale(const SimTK::State& s, const ScaleSet& scaleSet) override;

    const char* getWrapTypeName() const override;

    std::string getDimensionsString() const override;

    const double getSemiAxisLengthH() const;
    void setSemiAxisLengthH(double aSemiAxisLengthH);
    const double getSemiAxisLengthG() const;
    void setSemiAxisLengthG(double aSemiAxisLengthG);

protected:
    int wrapLine(const SimTK::State& s, SimTK::Vec3& aPoint1,
            SimTK::Vec3& aPoint2, const PathWrap& aPathWrap,
            WrapResult& aWrapResult, bool& aFlag) const override;

    int neglect_ellipse(SimTK::Vec3& aPoint1, SimTK::Vec3& aPoint2) const;

    /// Implement generateDecorations to draw geometry in visualizer
    void generateDecorations(bool fixed, const ModelDisplayHints& hints,
            const SimTK::State& state,
            SimTK::Array_<SimTK::DecorativeGeometry>& appendToThis)
            const override;

    void extendFinalizeFromProperties() override;

private:
    void constructProperties();

}; // END of class WrapViaEllipse

} // end of namespace OpenSim

#endif // OPENSIM_WRAP_HAMMER_ELLIPSE_H_
