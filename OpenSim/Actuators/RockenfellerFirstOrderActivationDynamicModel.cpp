/* -------------------------------------------------------------------------- *
 *            OpenSim:  RockenfellerFirstOrderActivationDynamicModel.cpp            *
 * -------------------------------------------------------------------------- *
*/
#include "RockenfellerFirstOrderActivationDynamicModel.h"

using namespace std;
using namespace OpenSim;
using namespace SimTK;

//==============================================================================
// CONSTRUCTION
//==============================================================================
RockenfellerFirstOrderActivationDynamicModel::RockenfellerFirstOrderActivationDynamicModel()
{
    setNull();
    constructProperties();
    setName("default_RockenfellerFirstOrderActivationDynamicModel");
}

RockenfellerFirstOrderActivationDynamicModel::
RockenfellerFirstOrderActivationDynamicModel(double time_constant_hatze,
                                            double nue,
                                            double roh_0,
                                            double gamma_C,
                                            double minimum_gamma,
                                            double minimum_activation,
                                            double optimal_fiber_length,
                                            const std::string& muscleName)
{
    setNull();
    constructProperties();

    std::string name = muscleName + "_activation";
    setName(name);
    
    set_time_constant_hatze(time_constant_hatze);
    set_nue(nue);
    set_roh_0(roh_0);
    set_gamma_C(gamma_C);
    set_minimum_gamma(minimum_gamma);
    set_minimum_activation(minimum_activation);
    set_optimal_fiber_length(optimal_fiber_length);
}

void RockenfellerFirstOrderActivationDynamicModel::setNull()
{
    setAuthors("Mike Spahr");
}

void RockenfellerFirstOrderActivationDynamicModel::constructProperties()
{
    constructProperty_time_constant_hatze(11.3);
    constructProperty_nue(3);
    constructProperty_roh_0(5.27);
    constructProperty_gamma_C(1.37);
    constructProperty_minimum_gamma(0.0);
    constructProperty_minimum_activation(0.005); // Hatze constant
    constructProperty_optimal_fiber_length(0.2);
}

//==============================================================================
// SERVICES
//==============================================================================
double RockenfellerFirstOrderActivationDynamicModel::
clampActivation(double activation) const
{
    return clamp(get_minimum_activation(), activation, 1.0);
}

double RockenfellerFirstOrderActivationDynamicModel::
calculateActivation(double gamma, double fiber_length) const
{
    double gamma_C = get_gamma_C();
    double rho_0 = get_roh_0();
    double nue = get_nue();
    double optimal_fiber_length = get_optimal_fiber_length();
    double min_act = get_minimum_activation();

    double rho = gamma_C * rho_0 * fiber_length / optimal_fiber_length;
    double rhogam = std::pow(gamma * rho, nue);
    return ((min_act + rhogam) / (1.0 + rhogam));
}

double RockenfellerFirstOrderActivationDynamicModel::
clampGamma(double gamma) const
{
    return clamp(get_minimum_gamma(), gamma, 1.0);
}

double RockenfellerFirstOrderActivationDynamicModel::
calcDerivative(double gamma, double excitation) const
{
    gamma = clamp(get_minimum_gamma(), gamma, 1.0);
    double hatze_constant = get_time_constant_hatze();
    return hatze_constant * (excitation - gamma);
}

//==============================================================================
// COMPONENT INTERFACE
//==============================================================================
void RockenfellerFirstOrderActivationDynamicModel::extendFinalizeFromProperties()
{
    Super::extendFinalizeFromProperties();

    std::string errorLocation = getName() + 
        " RockenfellerFirstOrderActivationDynamicModel::extendFinalizeFromProperties";

    // Ensure property values are within appropriate ranges.
    OPENSIM_THROW_IF_FRMOBJ(
        get_time_constant_hatze() < SimTK::SignificantReal,
        InvalidPropertyValue,
        getProperty_time_constant_hatze().getName(),
        "Time constant hatze must be greater than zero");
    OPENSIM_THROW_IF_FRMOBJ(
        get_nue() < SimTK::SignificantReal,
        InvalidPropertyValue,
        getProperty_nue().getName(),
        "Nue constant must be greater than zero");
    OPENSIM_THROW_IF_FRMOBJ(
        get_roh_0() < SimTK::SignificantReal,
        InvalidPropertyValue,
        getProperty_roh_0().getName(),
        "roh_0 constant must be greater than zero");
    OPENSIM_THROW_IF_FRMOBJ(
        get_gamma_C() < SimTK::SignificantReal,
        InvalidPropertyValue,
        getProperty_gamma_C().getName(),
        "gamma_C constant must be greater than zero");
    OPENSIM_THROW_IF_FRMOBJ(
        get_minimum_gamma() < 0 ||
        get_minimum_gamma() > 1.0-SimTK::SignificantReal,
        InvalidPropertyValue,
        getProperty_minimum_gamma().getName(),
        "Minimum calcium Concentration must be in the range [0, 1)");
}
