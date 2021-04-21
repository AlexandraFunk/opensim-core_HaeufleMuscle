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
    constructProperty_semi_axis_length_H(1.0); // TODO set realistic default value
    constructProperty_semi_axis_length_G(1.0); // TODO set realistic default value
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
