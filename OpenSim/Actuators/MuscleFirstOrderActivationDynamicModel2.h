#ifndef OPENSIM_MUSCLEFIRSTORDERACTIVATIONDYNAMICMODEL2_H_
#define OPENSIM_MUSCLEFIRSTORDERACTIVATIONDYNAMICMODEL2_H_
#include <OpenSim/Actuators/osimActuatorsDLL.h>
#include <OpenSim/Simulation/Model/ModelComponent.h>

namespace OpenSim {
class OSIMACTUATORS_API MuscleFirstOrderActivationDynamicModel2 : public ModelComponent{
OpenSim_DECLARE_CONCRETE_OBJECT(MuscleFirstOrderActivationDynamicModel2, ModelComponent);
public:

//==============================================================================
// PROPERTIES
//==============================================================================
    OpenSim_DECLARE_PROPERTY(minimum_activation, double,
        "Lower bound on activation (overridden when this is a subcomponent of a Muscle)");

    OpenSim_DECLARE_PROPERTY(activation_Hatze_time_constant, double,
            "Time constant, in 1/seconds. (overridden when this is a "
            "subcomponent of a Muscle)");
    OpenSim_DECLARE_PROPERTY(activation_exponent, double,
            "Hatze Coefficient (overridden when this is a subcomponent of a "
            "Muscle)");
    OpenSim_DECLARE_PROPERTY(activation_optimal_calcium_concentration_fraction,
            double, "roh_0 * gamma_C");
    OpenSim_DECLARE_PROPERTY(minimum_gamma, double,
            "Calcium concentration lower bound. (overridden when this is a "
            "subcomponent of a Muscle)");
    OpenSim_DECLARE_PROPERTY(optimal_fiber_length, double,
            "Optimal fiber length. (overridden when this is a subcomponent of "
            "a Muscle)");

//==============================================================================
// PUBLIC METHODS
//==============================================================================
    /** The default constructor creates an activation dynamic model with the
    default property values and assigns it a default name. **/
    MuscleFirstOrderActivationDynamicModel2();

    /** Creates an activation dynamic model using the provided parameters. */
    MuscleFirstOrderActivationDynamicModel2(
            double activation_Hatze_time_constant, double activation_exponent,
            double activation_optimal_calcium_concentration_fraction,
            double minimum_gamma, double optimal_fiber_length,
            double minActivation, const std::string& muscleName);

    /**
    @returns Activation clamped to the range [minActivation, 1.0].
    */
    double clampActivation(double activation) const;

    double calculateActivation(double gamma, double fiber_length) const;

    /**
     @returns Calcium concentration clamped to the range [minGamma, 1.0]
     */
    double clampGamma(double gamma) const;

    /** Calculates the time derivative of calcium concentration.. */
    double calcDerivative(double gamma, double excitation) const;

protected:
    // Component interface.
    void extendFinalizeFromProperties() override;

private:
    void setNull();
    void constructProperties();

};

}
#endif //OPENSIM_MUSCLEFIRSTORDERACTIVATIONDYNAMICMODEL2_H_
