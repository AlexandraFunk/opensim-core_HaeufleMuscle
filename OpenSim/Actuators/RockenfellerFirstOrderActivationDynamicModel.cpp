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
RockenfellerFirstOrderActivationDynamicModel(double activation_Hatze_time_constant,
                                            double activation_exponent,
                                            double activation_optimal_calcium_concentration_fraction,
                                            double minimum_gamma,
                                            double activation_minimum,
                                            double optimal_fiber_length,
                                            const std::string& muscleName)
{
    setNull();
    constructProperties();

    std::string name = muscleName + "_activation";
    setName(name);
    
    set_activation_Hatze_time_constant(activation_Hatze_time_constant);
    set_activation_exponent(activation_exponent);
    set_activation_optimal_calcium_concentration_fraction(
            activation_optimal_calcium_concentration_fraction);
    set_minimum_gamma(minimum_gamma);
    set_activation_minimum(activation_minimum);
    set_optimal_fiber_length(optimal_fiber_length);
}

void RockenfellerFirstOrderActivationDynamicModel::setNull()
{
    setAuthors("Mike Spahr");
}

void RockenfellerFirstOrderActivationDynamicModel::constructProperties()
{
    constructProperty_activation_Hatze_time_constant(11.3);
    constructProperty_activation_exponent(3);
    constructProperty_activation_optimal_calcium_concentration_fraction(5.27 * 1.37); // roh_0 * gamma_C
    constructProperty_minimum_gamma(0.0);
    constructProperty_activation_minimum(0.005); // Hatze constant
    constructProperty_optimal_fiber_length(0.2);
}

//==============================================================================
// SERVICES
//==============================================================================
double RockenfellerFirstOrderActivationDynamicModel::
clampActivation(double activation) const
{
    return clamp(get_activation_minimum(), activation, 1.0);
}

double RockenfellerFirstOrderActivationDynamicModel::
calculateActivation(double gamma, double fiber_length) const
{
    double pomega = get_activation_optimal_calcium_concentration_fraction();
    double activation_exponent = get_activation_exponent();
    double optimal_fiber_length = get_optimal_fiber_length();
    double min_act = get_activation_minimum();

    double rho = pomega * fiber_length / optimal_fiber_length;
    double rhogam = std::pow(gamma * rho, activation_exponent);
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
    double hatze_constant = get_activation_Hatze_time_constant();
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
        get_activation_Hatze_time_constant() < SimTK::SignificantReal,
        InvalidPropertyValue,
        getProperty_activation_Hatze_time_constant().getName(),
        "Time constant hatze must be greater than zero");
    OPENSIM_THROW_IF_FRMOBJ(
        get_activation_exponent() < SimTK::SignificantReal,
        InvalidPropertyValue,
        getProperty_activation_exponent().getName(),
        "activation_exponent constant must be greater than zero");
    OPENSIM_THROW_IF_FRMOBJ(
        get_activation_optimal_calcium_concentration_fraction() <
                SimTK::SignificantReal,
        InvalidPropertyValue,
        getProperty_activation_optimal_calcium_concentration_fraction()
                .getName(),
        "Roh_0 * Gamma_C constant must be greater than zero");
    OPENSIM_THROW_IF_FRMOBJ(
        get_minimum_gamma() < 0 ||
        get_minimum_gamma() > 1.0-SimTK::SignificantReal,
        InvalidPropertyValue,
        getProperty_minimum_gamma().getName(),
        "Minimum calcium Concentration must be in the range [0, 1)");
}
