/* -------------------------------------------------------------------------- *
 *                      OpenSim:  FatigableMuscle.cpp                         *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Copyright (c) 2005-2017 Stanford University and the Authors                *
 * Author(s): Peter Loan, Ajay Seth                                           *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied    *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

//=============================================================================
// INCLUDES
//=============================================================================
#include "RockenfellerMillardMuscle.h"

//=============================================================================
// STATICS
//=============================================================================
using namespace std;
using namespace OpenSim;

//=============================================================================
// CONSTRUCTOR(S) AND DESTRUCTOR
//=============================================================================
//_____________________________________________________________________________
/*
 * Default constructor
 */
RockenfellerMillardMuscle::RockenfellerMillardMuscle()
{
    constructProperties();
}

//_____________________________________________________________________________
/*
 * Constructor.
 */
RockenfellerMillardMuscle::RockenfellerMillardMuscle(const std::string &name, double maxIsometricForce, 
                       double optimalFiberLength, double tendonSlackLength,
                       double pennationAngle, double time_constant_hatze, 
                       double nue, double roh_0, double gamma_C) :
        Super(name, maxIsometricForce, optimalFiberLength, tendonSlackLength,
                pennationAngle)
{
    constructProperties();
    setTimeConstantHatze(time_constant_hatze);
    setNue(nue);
    setRoh_0(roh_0);
    setGamma_C(gamma_C);
}

//_____________________________________________________________________________
/*
 * Construct and initialize properties.
 * All properties are added to the property set. Once added, they can be
 * read in and written to files.
 */
void RockenfellerMillardMuscle::constructProperties()
{
    setAuthors("Maria Hammer, Mike Spahr");
    constructProperty_time_constant_hatze(11.3); // TODO realistic values
    constructProperty_nue(3); // TODO realistic values
    constructProperty_roh_0(5); // TODO realistic values
    constructProperty_gamma_C(1); // TODO realistic values
    //constructProperty_default_normalized_calcium_concentration(0.1); // TODO realistic values
}

// Define new states and their derivatives in the underlying system
void RockenfellerMillardMuscle::extendAddToSystem(SimTK::MultibodySystem& system) const
{
    // Allow Millard2012EquilibriumMuscle to add its states, before extending
    Super::extendAddToSystem(system);
    
    // // Now add the states
    // addStateVariable("normalized_calcium_concentration");
    // // and their corresponding derivatives
    // addCacheVariable("normalized_calcium_concentration_deriv", 0.0, SimTK::Stage::Dynamics);
    // std::cout << "Setting the states is done" << std::endl;
}

// void FatigableMuscle::extendAddToSystem(SimTK::MultibodySystem& system) const
// {
//     // Allow Millard2012EquilibriumMuscle to add its states, before extending
//     Super::extendAddToSystem(system);

//     // Now add the states necessary to implement the fatigable behavior
//     addStateVariable("target_activation");
//     addStateVariable("active_motor_units");
//     addStateVariable("fatigued_motor_units");
//     // and their corresponding derivatives
//     addCacheVariable("target_activation_deriv", 0.0, SimTK::Stage::Dynamics);
//     addCacheVariable("active_motor_units_deriv", 0.0, SimTK::Stage::Dynamics);
//     addCacheVariable("fatigued_motor_units_deriv", 0.0, SimTK::Stage::Dynamics);
// }

void RockenfellerMillardMuscle::extendInitStateFromProperties(SimTK::State& s) const
{
    Super::extendInitStateFromProperties(s);
    std::cout << "After extendInitStateFromProperties" << std::endl;
    //setNormalizedCalciumConcentration(s, getDefaultNormalizedCalciumConcentration());
    // double current_gamma = getNormalizedCalciumConcentration(s);
    // double current_lce = 0; //getFiberLength(s);
    // double nue = getNue();
    // double Kuh0 = getDefaultNormalizedCalciumConcentration();
    // std::cout << "Before extendInitStateFromProperties" << std::endl;

    // double rhogam = std::pow(current_gamma * rho(current_lce), nue);
    // double activation = (Kuh0 + rhogam) / (1.0 + rhogam);
    // if (!get_ignore_activation_dynamics()) {
    //     Super::setActivation(s, activation);
    // }
    //Super::setStateVariableValue(s, "activation", activation);
}



void RockenfellerMillardMuscle::extendSetPropertiesFromState(const SimTK::State& s)
{
    Super::extendSetPropertiesFromState(s);
    //setDefaultNormalizedCalciumConcentration(getNormalizedCalciumConcentration(s));
}

//--------------------------------------------------------------------------
// GET & SET Properties
//--------------------------------------------------------------------------
void RockenfellerMillardMuscle::setTimeConstantHatze(double aTimeConstantHatze)
{
    set_time_constant_hatze(aTimeConstantHatze);
}
void RockenfellerMillardMuscle::setNue(double aNue)
{
    set_nue(aNue);
}
void RockenfellerMillardMuscle::setRoh_0(double roh_0)
{
    set_roh_0(roh_0);
}
void RockenfellerMillardMuscle::setGamma_C(double gamma_C)
{
    set_gamma_C(gamma_C);
}

//void RockenfellerMillardMuscle::setDefaultNormalizedCalciumConcentration(double anormalized_calcium_concentration)
//{
//    set_default_normalized_calcium_concentration(anormalized_calcium_concentration);
//}

//--------------------------------------------------------------------------
// GET & SET States and their derivatives
//--------------------------------------------------------------------------
// double RockenfellerMillardMuscle::getNormalizedCalciumConcentration(const SimTK::State& s) const
// {  return getStateVariableValue(s, "normalized_calcium_concentration");}

// void RockenfellerMillardMuscle::setNormalizedCalciumConcentration(SimTK::State& s, 
//                                                            double normalizedCa) const
// {  setStateVariableValue(s, "normalized_calcium_concentration", normalizedCa);}

// double RockenfellerMillardMuscle::getNormalizedCalciumConcentrationDeriv(const SimTK::State& s) const
// {  return getStateVariableDerivativeValue(s, "normalized_calcium_concentration");}

// void RockenfellerMillardMuscle::setNormalizedCalciumConcentrationDeriv(const SimTK::State& s, 
//                                                 double normalizedCaDeriv) const
// {  setStateVariableDerivativeValue(s, "normalized_calcium_concentration", normalizedCaDeriv);}

double RockenfellerMillardMuscle::rho(double lce) const
    {
        double gamma_C = getGamma_C();
        double rho_0 = getRoh_0();
        double optimalFiberlength = getOptimalFiberLength();
        return (gamma_C * rho_0 * lce/ optimalFiberlength);
    }

//=============================================================================
// COMPUTATION
//=============================================================================
//_____________________________________________________________________________
/**
 * Compute the derivatives of the muscle states.
 *
 * @param s  system state
 */
void RockenfellerMillardMuscle::computeStateVariableDerivatives(const SimTK::State& s) const
{
    // Allow Super to assign any state derivative values for states it allocated
    Super::computeStateVariableDerivatives(s);
    
    int nd = getNumStateVariables();
    std::cout << "Number of State Variables are: " << nd << std::endl;
    
    
    // compute the normalized calcium concentration deriv
    double excitation = getExcitation(s);
    double mhatze = getTimeConstantHatze();
    double current_gamma = Super::getActivationModel().clampActivation(getStateVariableValue(s, "activation"));
    double normalizedCalciumConcentrationDeriv = mhatze * (excitation - current_gamma);

    Super::setStateVariableDerivativeValue(s, "activation", normalizedCalciumConcentrationDeriv);
    //setNormalizedCalciumConcentrationDeriv(s, normalizedCalciumConcentrationDeriv);
}


void RockenfellerMillardMuscle::
calcMuscleDynamicsInfo(const SimTK::State& s, MuscleDynamicsInfo& mdi) const
{
    try {
        // Get the quantities that we've already computed.
        const MuscleLengthInfo &mli = Super::getMuscleLengthInfo(s);
        const FiberVelocityInfo &mvi = Super::getFiberVelocityInfo(s);
        double fiberStateClamped = mvi.userDefinedVelocityExtras[0];

        // Get the properties of this muscle.
        double tendonSlackLen = Super::getTendonSlackLength();
        double optFiberLen    = Super::getOptimalFiberLength();
        double fiso           = Super::getMaxIsometricForce();
        //double penHeight      = penMdl.getParallelogramHeight();
        const TendonForceLengthCurve& fseCurve = Super::get_TendonForceLengthCurve();

        // Compute dynamic quantities.
        double a = SimTK::NaN;
        if(!get_ignore_activation_dynamics()) {
            a = getActivationModel().clampActivation(
                    getStateVariableValue(s, "activation"));
        } else {
            a = getActivationModel().clampActivation(getControl(s));
        }

        // bei rhogam mit a rechnen
        double current_gamma = a;
        double current_lce = Super::getFiberLength(s);
        double nue = getNue();
        double Kuh0 = Super::getDefaultActivation();
        double rhogam = std::pow(current_gamma * rho(current_lce), nue);
        //double activation = (Kuh0 + rhogam) / (1.0 + rhogam);
        a = (Kuh0 + rhogam) / (1.0 + rhogam);
            

        // Compute the stiffness of the muscle fiber.
        SimTK_ERRCHK_ALWAYS(mli.fiberLength > SimTK::SignificantReal,
            "calcMuscleDynamicsInfo",
            "The muscle fiber has a length of 0, causing a singularity");
        SimTK_ERRCHK_ALWAYS(mli.cosPennationAngle > SimTK::SignificantReal,
            "calcMuscleDynamicsInfo",
            "Pennation angle is 90 degrees, causing a singularity");

        double fm           = 0.0; //total fiber force
        double aFm          = 0.0; //active fiber force
        double p1Fm         = 0.0; //passive conservative fiber force
        double p2Fm         = 0.0; //passive non-conservative fiber force
        double pFm          = 0.0; //total passive fiber force
        double fmAT         = 0.0;
        double dFm_dlce     = 0.0;
        double dFmAT_dlceAT = 0.0;
        double dFt_dtl      = 0.0;
        double Ke           = 0.0;

        if(fiberStateClamped < 0.5) { //flag is set to 0.0 or 1.0
            SimTK::Vec4 fiberForceV;

            fiberForceV = calcFiberForce(fiso, a,
                                         mli.fiberActiveForceLengthMultiplier,
                                         mvi.fiberForceVelocityMultiplier,
                                         mli.fiberPassiveForceLengthMultiplier,
                                         mvi.normFiberVelocity);
            fm   = fiberForceV[0];
            aFm  = fiberForceV[1];
            p1Fm = fiberForceV[2];
            p2Fm = fiberForceV[3];
            pFm  = p1Fm + p2Fm;

            // Every configuration except the rigid tendon chooses a fiber
            // velocity that ensures that the fiber does not generate a
            // compressive force. Here, we must enforce that the fiber generates
            // only tensile forces by saturating the damping force generated by
            // the parallel element.
            if(get_ignore_tendon_compliance()) {
                if(fm < 0) {
                    fm   = 0.0;
                    p2Fm = -aFm - p1Fm;
                    pFm  = p1Fm + p2Fm;
                }
            }

            fmAT = fm * mli.cosPennationAngle;
            dFm_dlce = calcFiberStiffness(fiso, a,
                                          mvi.fiberForceVelocityMultiplier,
                                          mli.normFiberLength, optFiberLen);
            const double dFmAT_dlce =
                calc_DFiberForceAT_DFiberLength(fm, dFm_dlce, mli.fiberLength,
                                                mli.sinPennationAngle,
                                                mli.cosPennationAngle);
            dFmAT_dlceAT = calc_DFiberForceAT_DFiberLengthAT(dFmAT_dlce,
                mli.sinPennationAngle, mli.cosPennationAngle, mli.fiberLength);

            // Compute the stiffness of the tendon.
            if(!get_ignore_tendon_compliance()) {
                dFt_dtl = fseCurve.calcDerivative(mli.normTendonLength,1)
                          *(fiso/tendonSlackLen);

                // Compute the stiffness of the whole musculotendon actuator.
                if (abs(dFmAT_dlceAT*dFt_dtl) > 0.0
                    && abs(dFmAT_dlceAT+dFt_dtl) > SimTK::SignificantReal) {
                    Ke = (dFmAT_dlceAT*dFt_dtl)/(dFmAT_dlceAT+dFt_dtl);
                }
            } else {
                dFt_dtl = SimTK::Infinity;
                Ke = dFmAT_dlceAT;
            }
        }

        double fse = 0.0;
        if(!get_ignore_tendon_compliance()) {
            fse = fseCurve.calcValue(mli.normTendonLength);
        } else {
            fse = fmAT/fiso;
        }

        mdi.activation                = a;
        mdi.fiberForce                = fm;
        mdi.fiberForceAlongTendon     = fmAT;
        mdi.normFiberForce            = fm/fiso;
        mdi.activeFiberForce          = aFm;
        mdi.passiveFiberForce         = pFm;
        mdi.tendonForce               = fse*fiso;
        mdi.normTendonForce           = fse;
        mdi.fiberStiffness            = dFm_dlce;
        mdi.fiberStiffnessAlongTendon = dFmAT_dlceAT;
        mdi.tendonStiffness           = dFt_dtl;
        mdi.muscleStiffness           = Ke;

        // Verify that the derivative of system energy minus work is zero within
        // a reasonable numerical tolerance.
        //double dphidt       = mvi.pennationAngularVelocity;
        double dFibPEdt     = p1Fm*mvi.fiberVelocity; //only conservative part
                                                      //of passive fiber force
        double dTdnPEdt     = fse*fiso*mvi.tendonVelocity;
        double dFibWdt      = -(mdi.activeFiberForce+p2Fm)*mvi.fiberVelocity;
        double dmcldt       = getLengtheningSpeed(s);
        double dBoundaryWdt = mdi.tendonForce*dmcldt;

        //double dSysEdt = (dFibPEdt + dTdnPEdt) - dFibWdt - dBoundaryWdt;
        //double tol = sqrt(SimTK::Eps);

        // Populate the power entries.
        mdi.fiberActivePower  = dFibWdt;
        mdi.fiberPassivePower = -(dFibPEdt);
        mdi.tendonPower       = -dTdnPEdt;
        mdi.musclePower       = -dBoundaryWdt;

        // Store quantities unique to this Muscle: the passive conservative
        // (elastic) fiber force and the passive non-conservative (damping)
        // fiber force.
        SimTK::Vector dynExtras = SimTK::Vector(2);
        dynExtras[0] = p1Fm; //elastic
        dynExtras[1] = p2Fm; //damping
        mdi.userDefinedDynamicsExtras = dynExtras;

    } catch(const std::exception &x) {
        std::string msg = "Exception caught in Millard2012EquilibriumMuscle::"
                          "calcMuscleDynamicsInfo from " + getName() + "\n"
                          + x.what();
        cerr << msg << endl;
        throw OpenSim::Exception(msg);
    }
}


SimTK::Vec4 RockenfellerMillardMuscle::calcFiberForce(double fiso, double a,
        double fal, double fv, double fpe, double dlceN) const {
    double beta = Super::getFiberDamping();
    double fa = fiso * (a * fal * fv);
    double fp1 = fiso * fpe;
    double fp2 = fiso * beta * dlceN;
    double fm = fa + (fp1 + fp2);

    SimTK::Vec4 fiberF;
    fiberF[0] = fm;
    fiberF[1] = fa;
    fiberF[2] = fp1; // conservative passive force
    fiberF[3] = fp2; // non-conservative passive force
    return fiberF;
}

double RockenfellerMillardMuscle::calcFiberStiffness(
        double fiso, double a, double fv, double lceN, double optFibLen) const {
    const FiberForceLengthCurve& fpeCurve = Super::get_FiberForceLengthCurve();
    const ActiveForceLengthCurve& falCurve = Super::get_ActiveForceLengthCurve();
    double DlceN_Dlce = 1.0 / optFibLen;
    double Dfal_Dlce = falCurve.calcDerivative(lceN, 1) * DlceN_Dlce;
    double Dfpe_Dlce = fpeCurve.calcDerivative(lceN, 1) * DlceN_Dlce;

    // DFm_Dlce
    return fiso * (a * Dfal_Dlce * fv + Dfpe_Dlce);
} 

double RockenfellerMillardMuscle::calc_DFiberForceAT_DFiberLength(
                double fiberForce, double fiberStiffness, double lce,
                double sinPhi, double cosPhi) const 
{
    double Dphi_Dlce =
            Super::getPennationModel().calc_DPennationAngle_DfiberLength(lce);
    double Dcosphi_Dlce = -sinPhi * Dphi_Dlce;

    // The stiffness of the fiber along the direction of the tendon. For small
    // changes in length parallel to the fiber, this quantity is
    // D(FiberForceAlongTendon) / D(fiberLength)
    // dFmAT/dlce = d/dlce( fiso * (a *fal*fv + fpe + beta*dlceN)*cosPhi )
    return fiberStiffness * cosPhi + fiberForce * Dcosphi_Dlce;
}

double RockenfellerMillardMuscle::calc_DFiberForceAT_DFiberLengthAT(
        double dFmAT_d_lce, double sinPhi, double cosPhi, double lce) const 
{
    double dphi_d_lce =
            Super::getPennationModel().calc_DPennationAngle_DfiberLength(lce);

    // The change in length of the fiber length along the tendon.
    // lceAT = lce*cos(phi)
    double DlceAT_Dlce = cosPhi - lce * sinPhi * dphi_d_lce;

    // dFmAT/dlceAT = (dFmAT/dlce)*(1/(dlceAT/dlce))
    //              = dFmAT/dlceAT
    return dFmAT_d_lce * (1.0 / DlceAT_Dlce);
}
