#ifndef OPENSIM_ROCKENFELLERFIRSTORDERACTIVATIONDYNAMICMODEL_H_
#define OPENSIM_ROCKENFELLERFIRSTORDERACTIVATIONDYNAMICMODEL_H_
/* -------------------------------------------------------------------------- *
 *             OpenSim:  RockenfellerFirstOrderActivationDynamicModel.h       *
 * -------------------------------------------------------------------------- *
*/
#include "osimPluginDLL.h"
#include <OpenSim/Actuators/osimActuatorsDLL.h>
#include <OpenSim/Simulation/Model/ModelComponent.h>

namespace OpenSim {
class OSIMPLUGIN_API RockenfellerFirstOrderActivationDynamicModel
        : public ModelComponent {
    OpenSim_DECLARE_CONCRETE_OBJECT(RockenfellerFirstOrderActivationDynamicModel, ModelComponent);
public:

//==============================================================================
// PROPERTIES
//==============================================================================
    OpenSim_DECLARE_PROPERTY(activation_Hatze_time_constant, double,
        "Time constant, in 1/seconds. (overridden when this is a subcomponent of a Muscle)");
    OpenSim_DECLARE_PROPERTY(activation_exponent, double,
        "Hatze Coefficient (overridden when this is a subcomponent of a Muscle)" );
    OpenSim_DECLARE_PROPERTY(
            activation_optimal_calcium_concentration_fraction, double, "roh_0 * gamma_C");
    OpenSim_DECLARE_PROPERTY(minimum_gamma, double,
        "Calcium concentration lower bound. (overridden when this is a subcomponent of a Muscle)");
    OpenSim_DECLARE_PROPERTY(activation_minimum, double,
            "Activation lower bound (Hatze constant 0.005) equal to Kuh0");
    OpenSim_DECLARE_PROPERTY(optimal_fiber_length, double,
        "Optimal fiber length. (overridden when this is a subcomponent of a Muscle)");

//==============================================================================
// PUBLIC METHODS
//==============================================================================
    /** The default constructor creates an activation dynamic model with the
    default property values and assigns it a default name. **/
    RockenfellerFirstOrderActivationDynamicModel();

    /** Creates an activation dynamic model using the provided parameters. */
    RockenfellerFirstOrderActivationDynamicModel(
            double activation_Hatze_time_constant, double activation_exponent,
            double activation_optimal_calcium_concentration_fraction,
            double minimum_gamma, double activation_minimum,
            double optimal_fiber_length, const std::string& muscleName);

    /**
    @returns Activation clamped to the range [minActivation, 1.0].
    */
    double clampActivation(double activation) const;
    
    double calculateActivation(double gamma, double fiber_length) const;
    
    /**
     @returns Calcium concentration clamped to the range [minGamma, 1.0]
     */
    double clampGamma(double gamma) const;

    /** Calculates the time derivative of calcium concentration. */
    double calcDerivative(double gamma, double excitation) const;

protected:
    // Component interface.
    void extendFinalizeFromProperties() override;

private:
    void setNull();
    void constructProperties();

};

}
#endif //OPENSIM_ROCKENFELLERFIRSTORDERACTIVATIONDYNAMICMODEL_H_
