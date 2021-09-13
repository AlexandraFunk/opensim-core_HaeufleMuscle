#include "Millard2012Muscle.h"
#include <OpenSim/Simulation/Model/Model.h>

using namespace std;
using namespace OpenSim;
using namespace SimTK;

const string Millard2012Muscle::STATE_CALCIUM_CONCENTRATION =
        "calcium_concentration";
const string Millard2012Muscle::
    STATE_FIBER_LENGTH_NAME = "fiber_length";
const double MIN_NONZERO_DAMPING_COEFFICIENT = 0.001;

//==============================================================================
// PROPERTIES
//==============================================================================
void Millard2012Muscle::setNull()
{   setAuthors("Matthew Millard, Tom Uchida, Ajay Seth"); }

void Millard2012Muscle::constructProperties()
{
    constructProperty_fiber_damping(0.1); //damped model used by default
    constructProperty_default_calcium_concentration(0.0);
    constructProperty_default_fiber_length(getOptimalFiberLength());

    constructProperty_activation_Hatze_time_constant(11.3);
    constructProperty_activation_exponent(3.0);
    constructProperty_activation_optimal_calcium_concentration_fraction(
            5.27 * 1.37); // roh_0 * gamma_C
    constructProperty_minimum_activation(0.01);

    constructProperty_maximum_pennation_angle(acos(0.1));

    constructProperty_ActiveForceLengthCurve(ActiveForceLengthCurve());
    constructProperty_ForceVelocityCurve(ForceVelocityCurve());
    constructProperty_FiberForceLengthCurve(FiberForceLengthCurve());
    constructProperty_TendonForceLengthCurve(TendonForceLengthCurve());

    setMinControl(get_minimum_activation());
}

void Millard2012Muscle::extendFinalizeFromProperties()
{
    Super::extendFinalizeFromProperties();

    // Switch to undamped model if damping coefficient is small.
    if (get_fiber_damping() < MIN_NONZERO_DAMPING_COEFFICIENT) {
        set_fiber_damping(0.0);
        use_fiber_damping = false;
    } else {
        use_fiber_damping = true;
    }

    // Set the names of the muscle curves.
    const std::string& namePrefix = getName();

    ActiveForceLengthCurve& falCurve = upd_ActiveForceLengthCurve();
    falCurve.setName(namePrefix + "_ActiveForceLengthCurve");

    ForceVelocityCurve& fvCurve = upd_ForceVelocityCurve();
    fvCurve.setName(namePrefix + "_ForceVelocityCurve");

    FiberForceLengthCurve& fpeCurve = upd_FiberForceLengthCurve();
    fpeCurve.setName(namePrefix + "_FiberForceLengthCurve");

    TendonForceLengthCurve& fseCurve = upd_TendonForceLengthCurve();
    fseCurve.setName(namePrefix + "_TendonForceLengthCurve");

    // Include fiber damping in the model only if the damping coefficient is
    // larger than MIN_NONZERO_DAMPING_COEFFICIENT. This is done to ensure
    // we remain sufficiently far from the numerical singularity at beta=0.
    use_fiber_damping = (getFiberDamping() >= MIN_NONZERO_DAMPING_COEFFICIENT);

    // To initialize, we need to construct an *inverse* force-velocity curve
    // from the parameters of the force-velocity curve.
    double conSlopeAtVmax   = fvCurve.getConcentricSlopeAtVmax();
    double conSlopeNearVmax = fvCurve.getConcentricSlopeNearVmax();
    double isometricSlope   = fvCurve.getIsometricSlope();
    double eccSlopeAtVmax   = fvCurve.getEccentricSlopeAtVmax();
    double eccSlopeNearVmax = fvCurve.getEccentricSlopeNearVmax();
    double conCurviness     = fvCurve.getConcentricCurviness();
    double eccCurviness     = fvCurve.getEccentricCurviness();
    double eccForceMax      = fvCurve.getMaxEccentricVelocityForceMultiplier();

//Ensure that the slopes near vmax are postive, finite, and above zero
    OPENSIM_THROW_IF_FRMOBJ(conSlopeNearVmax < SimTK::SqrtEps,
            InvalidPropertyValue, 
            "ForceVelocityCurve:concentric_slope_near_vmax",
            "Slope near concentric vmax cannot be less than SimTK::SqrtEps"
            "(1.49e-8)");

    OPENSIM_THROW_IF_FRMOBJ(eccSlopeNearVmax < SimTK::SqrtEps,
            InvalidPropertyValue, 
            "ForceVelocityCurve:eccentric_slope_near_vmax",
            "Slope near eccentric vmax cannot be less than SimTK::SqrtEps"
            "(1.49e-8)");


    // A few parameters may need to be adjusted to avoid singularities 
    // if the classic Hill formulation is being used with an elastic tendon
    if(!get_ignore_tendon_compliance() && !use_fiber_damping) {
        // Compliant tendon with no damping.
        OPENSIM_THROW_IF_FRMOBJ(get_minimum_activation() < 0.01,
            InvalidPropertyValue, getProperty_minimum_activation().getName(),
            "Minimum activation cannot be less than 0.01 when using"
            "the classic Hill model with an elastic tendon due to a"
            "singularity in the state derivative for an activation of 0.");

        OPENSIM_THROW_IF_FRMOBJ(getMinControl() < get_minimum_activation(),
            InvalidPropertyValue, getProperty_min_control().getName(),
            "Minimum control cannot be less than minimum activation");

        if (falCurve.getMinValue() < 0.1){
            log_info("'{}' Parameter update for classic Hill model: "
                "ActiveForceLengthCurve parameter '{}' was {} but is now {}.", 
                getName(),"minimum_value", falCurve.getMinValue(), 0.1);
            falCurve.setMinValue(0.1);
        }
        

        if (conSlopeAtVmax < conSlopeNearVmax*0.5){
            log_info("'{}' Parameter update for classic Hill model: "
                "ForceVelocityCurve parameter '{}' was {} but is now {}.", 
                getName(),"eccentric_slope_near_vmax", conSlopeAtVmax, 
                conSlopeNearVmax*0.5);            
            conSlopeAtVmax  = conSlopeNearVmax*0.5;            
        }
        if(conSlopeNearVmax < 0.05){               
            log_info("'{}': Warning slow simulation: classic Hill model"
                "is being used with a '{}' of {}. Use the damped model, or "
                "increase to {} for faster simulations.", getName(),
                "concentric_slope_near_vmax",conSlopeNearVmax, 0.05);                 
        }


        if (eccSlopeAtVmax < eccSlopeNearVmax*0.5){
            log_info("'{}' Parameter update for classic Hill model: "
                "ForceVelocityCurve parameter '{}' was {} but is now {}.", 
                getName(),"eccentric_slope_near_vmax", eccSlopeAtVmax, 
                eccSlopeNearVmax*0.5);
            eccSlopeAtVmax   = eccSlopeNearVmax*0.5;

        }
        if(eccSlopeNearVmax < 0.05){
            log_info("'{}': Warning slow simulation: classic Hill model"
                "is being used with a '{}' of {}. Use the damped model, or "
                "increase to {} for faster simulations.", getName(),
                "eccentric_slope_near_vmax",eccSlopeNearVmax, 0.05);                
        }        

        fvCurve.setCurveShape(
                conSlopeAtVmax, conSlopeNearVmax, isometricSlope,
                eccSlopeAtVmax, eccSlopeNearVmax, eccForceMax);  


    } else { 
        const double min_activation = get_minimum_activation();
        const double min_activation_clamped = clamp(0, min_activation, 1);
        if (min_activation > 0. && 
            std::abs(min_activation_clamped -min_activation) > SimTK::Eps) {
            log_info("'{}': Parameter update for the damped-model: "
                "minimum_activation was {} but is now {}",
                    getName(), min_activation, 
                   min_activation_clamped);

            set_minimum_activation(min_activation_clamped);
        }
        if(falCurve.getMinValue() > 0.0){
            log_info("'{}' Parameter update for the damped-model: "
                "ActiveForceLengthCurve minimum value was {} but is now {}.",
                 getName(), falCurve.getMinValue(), 0.);
            falCurve.setMinValue(0.0);
        }

        if(conSlopeAtVmax > 0.){
            log_info("'{}' Parameter update for the damped-model:"
                " ForceVelocityCurve  '{}' was {} but is now {}.", 
                getName(),"concentric_slope_at_vmax", conSlopeAtVmax, 0.);                        
            conSlopeAtVmax  = 0.;
        }
        if(eccSlopeAtVmax > 0.){            
            log_info("'{}' Parameter update for the damped-model: "
                "ForceVelocityCurve '{}' was {} but is now {}.", 
                getName(),"eccentric_slope_at_vmax", eccSlopeAtVmax, 0.);            
            eccSlopeAtVmax   = 0.;
        }

        fvCurve.setCurveShape(
            conSlopeAtVmax, conSlopeNearVmax, isometricSlope,
            eccSlopeAtVmax, eccSlopeNearVmax, eccForceMax);

        OPENSIM_THROW_IF_FRMOBJ(get_minimum_activation() < 0.0,
            InvalidPropertyValue, getProperty_minimum_activation().getName(),
            "Minimum activation cannot be less than zero");
        
        //singularity-free model still cannot have excitations below 0
        OPENSIM_THROW_IF_FRMOBJ(getMinControl() < get_minimum_activation(),
            InvalidPropertyValue, getProperty_min_control().getName(),
            "Minimum control cannot be less than minimum activation");
    }


    //The ForceVelocityInverseCurve is used with both the damped and the 
    // classic Hill model but in different ways:
    //
    // The damped model uses the fvInvCurve to provide an initial fiber 
    // velocity prior to numerically solving the root using Eqn. 8 from 
    // Millard et al.
    //
    // The classic Hill formulation uses the fvInvCurve to directly solve for 
    // fiber velocity using Eqn. 6 from Millard et al.
    //
    //While the fvCurve and fvInvCurve are inverses of each other between
    //concentric-velocity-near-vmax to eccentric-velocity-near-vmax these
    //curves differ outside of this region. Why? When an elastic-tendon model
    //is used with the classic Hill formulation Eqn. 6 of Millard et al. is 
    //used to evaluate the state derivative of fiber length. As noted in the 
    //paper Eqn. 6 has several singularites one of which is related to the 
    //force-velocity-inverse curve: simply, the force-velocity-inverse curve 
    //must be invertible.
    //
    //In this case this is equivalent to saying that the force-velocity-inverse 
    //curve must have finite endslopes. To avoid this singularity we set the 
    //end slopes of the force-velocity-inverse curve to be a fraction of the 
    //`near' slopes used by the force-velocity curve: this value is guaranteed
    //to be finite, and is also guaranteed not to break the monotonicity of the 
    //curve. Also don't forget: the ForceVelocityInverseCurve constructor takes 
    //in the parameters of the force-velocity-curve that it inverts. So if you 
    //pass in an endslope of 0 the inverse curve will have an endslope of 
    //1/0 = infinity! Instead we make sure that we pass in a finite value for 
    //the end slopes so that the inverse curve has finite endslopes. 
    //
    //The fvCurve is used exclusively by the damped model. As such the fvCurve
    //does not need to be invertible since the damped model's state derivative
    //(Eqn. 8 of Millard et al.) has no singularities. This means that the user
    //can set the endslopes to zero which is particularly useful if simulations
    //of very rapid contractions are being performed. This also means that the
    //fvInvCurve is not the inverse of fvCurve for concentric velocities faster 
    //than concentric-velocity-near-vmax and eccentric velocities faster than 
    //eccentric-velocity-near-vmax: the end slopes differ.
    //
    //Millard M, Uchida T, Seth A, Delp SL. Flexing computational muscle: 
    //modeling and simulation of musculotendon dynamics. Journal of 
    //biomechanical engineering. 2013 Feb 1;135(2).
    double eccSlopeAtVmaxFvInv  = eccSlopeNearVmax*0.5;
    double conSlopeAtVmaxFvInv  = conSlopeNearVmax*0.5;

    fvInvCurve = ForceVelocityInverseCurve(conSlopeAtVmaxFvInv, conSlopeNearVmax,
                                           isometricSlope, eccSlopeAtVmaxFvInv,
                                           eccSlopeNearVmax, eccForceMax,
                                           conCurviness, eccCurviness);

    // Ensure all muscle curves are up-to-date.
    falCurve.ensureCurveUpToDate();
    fvCurve.ensureCurveUpToDate();
    fvInvCurve.ensureCurveUpToDate();
    fpeCurve.ensureCurveUpToDate();
    fseCurve.ensureCurveUpToDate();

    // Propagate properties down to pennation model subcomponent. If any of the
    // new property values are invalid, restore the subcomponent's current
    // property values (to avoid throwing again when the subcomponent's
    // extendFinalizeFromProperties() method is called directly) and then
    // re-throw the exception thrown by the subcomponent.
    auto& penMdl =
        updMemberSubcomponent<MuscleFixedWidthPennationModel>(penMdlIdx);
    MuscleFixedWidthPennationModel penMdlCopy(penMdl);
    penMdl.set_optimal_fiber_length(getOptimalFiberLength());
    penMdl.set_pennation_angle_at_optimal(getPennationAngleAtOptimalFiberLength());
    penMdl.set_maximum_pennation_angle(get_maximum_pennation_angle());
    try {
        penMdl.finalizeFromProperties();
    } catch (const InvalidPropertyValue&) {
        penMdl = penMdlCopy;
        throw;
    }

    // Propagate properties down to activation dynamics model subcomponent.
    // Handle invalid properties as above for pennation model.
    if (!get_ignore_activation_dynamics()) {
        auto& actMdl =
            updMemberSubcomponent<MuscleFirstOrderActivationDynamicModel2>(actMdlIdx);
        MuscleFirstOrderActivationDynamicModel2 actMdlCopy(actMdl);
        actMdl.set_optimal_fiber_length(get_optimal_fiber_length());
        actMdl.set_minimum_gamma(get_default_calcium_concentration());
        actMdl.set_activation_Hatze_time_constant(
                get_activation_Hatze_time_constant());
        actMdl.set_activation_exponent(get_activation_exponent());
        actMdl.set_activation_optimal_calcium_concentration_fraction(
                get_activation_optimal_calcium_concentration_fraction());
        actMdl.set_minimum_activation(get_minimum_activation());
        try {
            actMdl.finalizeFromProperties();
        } catch (const InvalidPropertyValue&) {
            actMdl = actMdlCopy;
            throw;
        }
    }

    // Compute and store values that are used for clamping the fiber length.
    const double minActiveFiberLength = falCurve.getMinActiveFiberLength()
                                        * getOptimalFiberLength();
    const double minPennatedFiberLength = penMdl.getMinimumFiberLength();
    m_minimumFiberLength = max(SimTK::SignificantReal,
        max(minActiveFiberLength, minPennatedFiberLength));

    const double phi = penMdl.calcPennationAngle(m_minimumFiberLength);
    m_minimumFiberLengthAlongTendon =
        penMdl.calcFiberLengthAlongTendon(m_minimumFiberLength, cos(phi));
}

//==============================================================================
// CONSTRUCTORS
//==============================================================================
Millard2012Muscle::Millard2012Muscle()
{
    setNull();
    constructProperties();
}

Millard2012Muscle::Millard2012Muscle(
const std::string &aName, double aMaxIsometricForce, double aOptimalFiberLength,
double aTendonSlackLength, double aPennationAngle)
{
    setNull();
    constructProperties();

    setName(aName);
    setMaxIsometricForce(aMaxIsometricForce);
    setOptimalFiberLength(aOptimalFiberLength);
    setTendonSlackLength(aTendonSlackLength);
    setPennationAngleAtOptimalFiberLength(aPennationAngle);
}

//==============================================================================
// GET METHODS
//==============================================================================
double Millard2012Muscle::getFiberDamping() const
{   return get_fiber_damping(); }
double Millard2012Muscle::getDefaultCalciumConcentration() const {
    return get_default_calcium_concentration();
}
double Millard2012Muscle::getDefaultFiberLength() const
{   return get_default_fiber_length(); }


const ActiveForceLengthCurve& Millard2012Muscle::
getActiveForceLengthCurve() const
{   return get_ActiveForceLengthCurve(); }

const ForceVelocityCurve& Millard2012Muscle::
getForceVelocityCurve() const
{   return get_ForceVelocityCurve(); }

const FiberForceLengthCurve& Millard2012Muscle::
getFiberForceLengthCurve() const
{   return get_FiberForceLengthCurve(); }

const TendonForceLengthCurve& Millard2012Muscle::
getTendonForceLengthCurve() const
{   return get_TendonForceLengthCurve(); }

const MuscleFixedWidthPennationModel& Millard2012Muscle::
getPennationModel() const
{ return getMemberSubcomponent<MuscleFixedWidthPennationModel>(penMdlIdx); }

const MuscleFirstOrderActivationDynamicModel2& Millard2012Muscle::
getActivationModel() const
{ return getMemberSubcomponent<MuscleFirstOrderActivationDynamicModel2>(actMdlIdx); }

double Millard2012Muscle::getMinimumFiberLength() const
{   return m_minimumFiberLength; }
double Millard2012Muscle::getMinimumFiberLengthAlongTendon() const
{   return m_minimumFiberLengthAlongTendon; }

double Millard2012Muscle::
getTendonForceMultiplier(SimTK::State& s) const
{   return getMuscleDynamicsInfo(s).normTendonForce; }

double Millard2012Muscle::
getFiberStiffnessAlongTendon(const SimTK::State& s) const
{ return getMuscleDynamicsInfo(s).fiberStiffnessAlongTendon; }

double Millard2012Muscle::
getFiberVelocity(const SimTK::State& s) const
{   return getFiberVelocityInfo(s).fiberVelocity; }

double Millard2012Muscle::getCalciumConcentration(const SimTK::State& s) const {
    return getStateVariableValue(s, STATE_CALCIUM_CONCENTRATION);
}

double Millard2012Muscle::getCalciumConcentrationDerivative(
        const SimTK::State& s) const {
    if (get_ignore_activation_dynamics()) { return 0.0; }

    return getActivationModel().calcDerivative(
            getCalciumConcentration(s), getExcitation(s));
}

double Millard2012Muscle::
getPassiveFiberElasticForce(const SimTK::State& s) const
{
    return getMuscleDynamicsInfo(s).userDefinedDynamicsExtras[0];
}

double Millard2012Muscle::
getPassiveFiberElasticForceAlongTendon(const SimTK::State& s) const
{
    return getMuscleDynamicsInfo(s).userDefinedDynamicsExtras[0] *
           getMuscleLengthInfo(s).cosPennationAngle;
}

double Millard2012Muscle::
getPassiveFiberDampingForce(const SimTK::State& s) const
{
    return getMuscleDynamicsInfo(s).userDefinedDynamicsExtras[1];
}

double Millard2012Muscle::
getPassiveFiberDampingForceAlongTendon(const SimTK::State& s) const
{
    return getMuscleDynamicsInfo(s).userDefinedDynamicsExtras[1] *
           getMuscleLengthInfo(s).cosPennationAngle;
}


//==============================================================================
// SET METHODS
//==============================================================================
void Millard2012Muscle::
setMuscleConfiguration(bool ignoreTendonCompliance,
                       bool ignoreActivationDynamics,
                       double dampingCoefficient)
{
    set_ignore_tendon_compliance(ignoreTendonCompliance);
    set_ignore_activation_dynamics(ignoreActivationDynamics);
    setFiberDamping(dampingCoefficient);
    finalizeFromProperties();
}

void Millard2012Muscle::setFiberDamping(double dampingCoefficient)
{   set_fiber_damping(dampingCoefficient); }

void Millard2012Muscle::setDefaultCalciumConcentration(
        double calciumConcentration) {
    set_default_calcium_concentration(calciumConcentration);
}

void Millard2012Muscle::
setActivation(SimTK::State& s, double activation) const
{
    if (get_ignore_activation_dynamics()) {
        SimTK::Vector& controls(_model->updControls(s));
        setControls(SimTK::Vector(1, activation), controls);
        _model->setControls(s, controls);
    } else {
        setStateVariableValue(s, STATE_CALCIUM_CONCENTRATION,
                getActivationModel().clampGamma(activation));
    }
    markCacheVariableInvalid(s, _velInfoCV);
    markCacheVariableInvalid(s, _dynamicsInfoCV);
}

void Millard2012Muscle::setDefaultFiberLength(double fiberLength)
{   set_default_fiber_length(fiberLength); }

void Millard2012Muscle::
setMinimumActivation(double minimumActivation)
{   set_minimum_activation(minimumActivation); }

void Millard2012Muscle::setActiveForceLengthCurve(
ActiveForceLengthCurve& aActiveForceLengthCurve)
{   set_ActiveForceLengthCurve(aActiveForceLengthCurve); }

void Millard2012Muscle::setForceVelocityCurve(
ForceVelocityCurve& aForceVelocityCurve)
{   set_ForceVelocityCurve(aForceVelocityCurve); }

void Millard2012Muscle::setFiberForceLengthCurve(
FiberForceLengthCurve& aFiberForceLengthCurve)
{   set_FiberForceLengthCurve(aFiberForceLengthCurve); }

void Millard2012Muscle::setTendonForceLengthCurve(
TendonForceLengthCurve& aTendonForceLengthCurve)
{   set_TendonForceLengthCurve(aTendonForceLengthCurve); }

void Millard2012Muscle::
setFiberLength(SimTK::State& s, double fiberLength) const
{
    if (!get_ignore_tendon_compliance()) {
        setStateVariableValue(s, STATE_FIBER_LENGTH_NAME,
                              clampFiberLength(fiberLength));
        markCacheVariableInvalid(s, _lengthInfoCV);
        markCacheVariableInvalid(s, _velInfoCV);
        markCacheVariableInvalid(s, _dynamicsInfoCV);
    }
}

//==============================================================================
// MUSCLE.H INTERFACE
//==============================================================================
double Millard2012Muscle::
computeActuation(const SimTK::State& s) const
{
    const MuscleDynamicsInfo& mdi = getMuscleDynamicsInfo(s);
    setActuation(s, mdi.tendonForce);
    return mdi.tendonForce;
}


void Millard2012Muscle::
computeFiberEquilibrium(SimTK::State& s, bool solveForVelocity) const
{
    if(get_ignore_tendon_compliance()) {                    // rigid tendon
        return;
    }

    // Elastic tendon initialization routine.

    // Initialize activation and fiber length provided by the State s
    _model->getMultibodySystem().realize(s, SimTK::Stage::Velocity);

    // Compute the fiber length where the fiber and tendon are in static
    // equilibrium. Fiber and tendon velocity are set to zero.

    // tol is the desired tolerance in Newtons.
    const double tol = max(1e-8*getMaxIsometricForce(), SimTK::SignificantReal*10);

    int maxIter = 200;
    double pathLength = getLength(s);
    double pathSpeed = solveForVelocity ? getLengtheningSpeed(s) : 0;
    double activation = getActivation(s);

    try {
        std::pair<StatusFromEstimateMuscleFiberState,
                  ValuesFromEstimateMuscleFiberState> result =
            estimateMuscleFiberState(activation, pathLength, pathSpeed,
                tol, maxIter, solveForVelocity);

        switch(result.first) {

        case StatusFromEstimateMuscleFiberState::Success_Converged:
            setActuation(s, result.second["tendon_force"]);
            setFiberLength(s, result.second["fiber_length"]);
            break;

        case StatusFromEstimateMuscleFiberState::Warning_FiberAtLowerBound:
            log_warn("Millard2012Muscle static solution: '{}' is "
                   "at its minimum fiber length of {}.",
                   getName(), result.second["fiber_length"]);
            setActuation(s, result.second["tendon_force"]);
            setFiberLength(s, result.second["fiber_length"]);
            break;

        case StatusFromEstimateMuscleFiberState::Failure_MaxIterationsReached:
            // Report internal variables and throw exception.
            std::ostringstream ss;
            ss << "\n  Solution error " << abs(result.second["solution_error"])
               << " exceeds tolerance of " << tol << "\n"
               << "  Newton iterations reached limit of " << maxIter << "\n"
               << "  Activation is " << activation << "\n"
               << "  Fiber length is " << result.second["fiber_length"] << "\n";
            OPENSIM_THROW_FRMOBJ(MuscleCannotEquilibrate, ss.str());
            break;
        }

    } catch (const std::exception& x) {
        OPENSIM_THROW_FRMOBJ(MuscleCannotEquilibrate,
            "Internal exception encountered.\n" + std::string{x.what()});
    }
}

//==============================================================================
// SCALING
//==============================================================================
void Millard2012Muscle::
extendPostScale(const SimTK::State& s, const ScaleSet& scaleSet)
{
    Super::extendPostScale(s, scaleSet);

    GeometryPath& path = upd_GeometryPath();
    if (path.getPreScaleLength(s) > 0.0)
    {
        double scaleFactor = path.getLength(s) / path.getPreScaleLength(s);
        upd_optimal_fiber_length() *= scaleFactor;
        upd_tendon_slack_length() *= scaleFactor;

        // Clear the pre-scale length that was stored in the GeometryPath.
        path.setPreScaleLength(s, 0.0);
    }
}

//==============================================================================
// MUSCLE INTERFACE REQUIREMENTS -- MUSCLE LENGTH INFO
//==============================================================================
void Millard2012Muscle::calcMuscleLengthInfo(const SimTK::State& s,
    MuscleLengthInfo& mli) const
{
    // Get musculotendon actuator properties.
    //double maxIsoForce    = getMaxIsometricForce();
    double optFiberLength = getOptimalFiberLength();
    double tendonSlackLen = getTendonSlackLength();

    try {
        // Get muscle-specific properties.
        const ActiveForceLengthCurve& falCurve = get_ActiveForceLengthCurve();
        const FiberForceLengthCurve&  fpeCurve = get_FiberForceLengthCurve();
        //const TendonForceLengthCurve& fseCurve = get_TendonForceLengthCurve();

        if(get_ignore_tendon_compliance()) {                //rigid tendon
            mli.fiberLength = clampFiberLength(
                               getPennationModel().calcFiberLength(getLength(s),
                               tendonSlackLen));
        } else {                                            // elastic tendon
            mli.fiberLength = clampFiberLength(
                                getStateVariableValue(s, STATE_FIBER_LENGTH_NAME));
        }

        mli.normFiberLength   = mli.fiberLength / optFiberLength;
        mli.pennationAngle    = getPennationModel().
                                    calcPennationAngle(mli.fiberLength);
        mli.cosPennationAngle = cos(mli.pennationAngle);
        mli.sinPennationAngle = sin(mli.pennationAngle);
        mli.fiberLengthAlongTendon = mli.fiberLength * mli.cosPennationAngle;

        // Necessary even for the rigid tendon, as it might have gone slack.
        mli.tendonLength      = getPennationModel().
                                    calcTendonLength(mli.cosPennationAngle,
                                                     mli.fiberLength,
                                                     getLength(s));
        mli.normTendonLength  = mli.tendonLength / tendonSlackLen;
        mli.tendonStrain      = mli.normTendonLength - 1.0;

        mli.fiberPassiveForceLengthMultiplier =
            fpeCurve.calcValue(mli.normFiberLength);
        mli.fiberActiveForceLengthMultiplier =
            falCurve.calcValue(mli.normFiberLength);

    } catch(const std::exception &x) {
        std::string msg = "Exception caught in Millard2012Muscle::"
                          "calcMuscleLengthInfo from " + getName() + "\n"
                          + x.what();
        throw OpenSim::Exception(msg);
    }
}


//==============================================================================
// MUSCLE INTERFACE REQUIREMENTS -- MUSCLE POTENTIAL ENERGY INFO
//==============================================================================
void Millard2012Muscle::
    calcMusclePotentialEnergyInfo(const SimTK::State& s,
        MusclePotentialEnergyInfo& mpei) const
{
    // Get musculotendon actuator properties.
    double maxIsoForce    = getMaxIsometricForce();
    double optFiberLength = getOptimalFiberLength();
    double tendonSlackLen = getTendonSlackLength();

    try {
        // Get the quantities that we've already computed.
        const MuscleLengthInfo &mli = getMuscleLengthInfo(s);

        // Get muscle-specific properties.
        //const ActiveForceLengthCurve& falCurve = get_ActiveForceLengthCurve();
        const FiberForceLengthCurve&  fpeCurve = get_FiberForceLengthCurve();
        const TendonForceLengthCurve& fseCurve = get_TendonForceLengthCurve();


        // Note that the curves return normalized area.
        double fiberStrainAtFiso     = fpeCurve.getStrainAtOneNormForce()
                                       - fpeCurve.getStrainAtZeroForce();
        double fiberStretchAtFiso    = fiberStrainAtFiso * optFiberLength;
        double fiberPotentialScaling = (fiberStretchAtFiso*maxIsoForce)
                                       / fiberStrainAtFiso;

        mpei.fiberPotentialEnergy = fpeCurve.calcIntegral(mli.normFiberLength)
                                   * fiberPotentialScaling;

        mpei.tendonPotentialEnergy = 0;
        if(!get_ignore_tendon_compliance()) {               // elastic tendon
            double tendonStrainAtFiso  = fseCurve.getStrainAtOneNormForce();
            double tendonStretchAtFiso = tendonStrainAtFiso*tendonSlackLen;
            double tendonPotentialScaling = (tendonStretchAtFiso*maxIsoForce)
                                            / tendonStrainAtFiso;
            mpei.tendonPotentialEnergy =
                fseCurve.calcIntegral(mli.normTendonLength)
                * tendonPotentialScaling;
        }

        mpei.musclePotentialEnergy = mpei.fiberPotentialEnergy +
                                     mpei.tendonPotentialEnergy;

    } catch(const std::exception &x) {
        std::string msg = "Exception caught in Millard2012Muscle::"
                          "calcMusclePotentialEnergyInfo from " + getName() + "\n"
                          + x.what();
        throw OpenSim::Exception(msg);
    }
}


//==============================================================================
// MUSCLE INTERFACE REQUIREMENTS -- FIBER VELOCITY INFO
//==============================================================================
void Millard2012Muscle::
calcFiberVelocityInfo(const SimTK::State& s, FiberVelocityInfo& fvi) const
{
    try {
        // Get the quantities that we've already computed.
        const MuscleLengthInfo &mli = getMuscleLengthInfo(s);

        // Get the static properties of this muscle.
        double dlenMcl   = getLengtheningSpeed(s);
        double optFibLen = getOptimalFiberLength();

        //======================================================================
        // Compute fv by inverting the force-velocity relationship in the
        // equilibrium equations.
        //======================================================================
        double dlce  = SimTK::NaN;
        double dlceN = SimTK::NaN;
        double fv    = SimTK::NaN;

        // Calculate fiber velocity.
        if(get_ignore_tendon_compliance()) {

            // Rigid tendon.

            if(mli.tendonLength < getTendonSlackLength()
                                  - SimTK::SignificantReal) {
                // The tendon is buckling, so fiber velocity is zero.
                dlce  = 0.0;
                dlceN = 0.0;
                fv    = 1.0;
            } else {
                dlce = getPennationModel().
                         calcFiberVelocity(mli.cosPennationAngle, dlenMcl, 0.0);
                dlceN = dlce/(optFibLen*getMaxContractionVelocity());
                fv = get_ForceVelocityCurve().calcValue(dlceN);
            }

        } else if(!get_ignore_tendon_compliance() && !use_fiber_damping) {

            // Elastic tendon, no damping.

            double a = SimTK::NaN;
            if(!get_ignore_activation_dynamics()) {
                a = getActivationModel().clampActivation(
                        getActivationModel().calculateActivation(
                                getStateVariableValue(
                                        s, STATE_CALCIUM_CONCENTRATION),
                                mli.fiberLength));
            } else {
                a = getActivationModel().clampActivation(getControl(s));
            }

            const TendonForceLengthCurve& fseCurve =
                get_TendonForceLengthCurve();
            double fse = fseCurve.calcValue(mli.normTendonLength);

            SimTK_ERRCHK_ALWAYS(mli.cosPennationAngle > SimTK::SignificantReal,
                "calcFiberVelocityInfo",
                "%s: Pennation angle is 90 degrees, causing a singularity");
            SimTK_ERRCHK_ALWAYS(a > SimTK::SignificantReal,
                "calcFiberVelocityInfo",
                "%s: Activation is 0, causing a singularity");
            SimTK_ERRCHK_ALWAYS(mli.fiberActiveForceLengthMultiplier >
                                SimTK::SignificantReal,
                "calcFiberVelocityInfo",
                "%s: Active-force-length factor is 0, causing a singularity");

            fv = calcFv(a, mli.fiberActiveForceLengthMultiplier,
                        mli.fiberPassiveForceLengthMultiplier, fse,
                        mli.cosPennationAngle);

            // Evaluate the inverse force-velocity curve.
            dlceN = fvInvCurve.calcValue(fv);
            dlce  = dlceN*getMaxContractionVelocity()*optFibLen;

        } else {

            // Elastic tendon, with damping.

            double a = SimTK::NaN;
            if(!get_ignore_activation_dynamics()) {
                a = getActivationModel().clampActivation(
                        getActivationModel().calculateActivation(
                                getStateVariableValue(
                                        s, STATE_CALCIUM_CONCENTRATION),
                                mli.fiberLength));
            } else {
                a = getActivationModel().clampActivation(getControl(s));
            }

            const TendonForceLengthCurve& fseCurve =
                get_TendonForceLengthCurve();
            double fse = fseCurve.calcValue(mli.normTendonLength);

            // Newton solve for fiber velocity.
            fv = 1.0;
            dlce = -1;
            dlceN = -1;
            double beta = get_fiber_damping();

            SimTK_ERRCHK_ALWAYS(beta > SimTK::SignificantReal,
                "calcFiberVelocityInfo",
                "Fiber damping coefficient must be greater than 0.");

            SimTK::Vec3 fiberVelocityV = calcDampedNormFiberVelocity(
                getMaxIsometricForce(), a, mli.fiberActiveForceLengthMultiplier,
                mli.fiberPassiveForceLengthMultiplier, fse, beta,
                mli.cosPennationAngle);

            // If the Newton method converged, update the fiber velocity.
            if(fiberVelocityV[2] > 0.5) { //flag is set to 0.0 or 1.0
                dlceN = fiberVelocityV[0];
                dlce  = dlceN*getOptimalFiberLength()
                        *getMaxContractionVelocity();
                fv = get_ForceVelocityCurve().calcValue(dlceN);
            } else {
                // Throw an exception here because there is no point integrating
                // a muscle velocity that is invalid (it will end up producing
                // invalid fiber lengths and will ultimately cause numerical
                // problems). The idea is to produce an exception and catch this
                // early before it can cause more damage.
                throw (OpenSim::Exception(getName() +
                       " Fiber velocity Newton method did not converge"));
            }
        }

        // Compute the other velocity-related components.
        double dphidt = getPennationModel().calcPennationAngularVelocity(
            tan(mli.pennationAngle), mli.fiberLength, dlce);
        double dlceAT = getPennationModel().calcFiberVelocityAlongTendon(
            mli.fiberLength, dlce, mli.sinPennationAngle, mli.cosPennationAngle,
            dphidt);
        double dmcldt = getLengtheningSpeed(s);
        double dtl = 0;

        if(!get_ignore_tendon_compliance()) {
            dtl = getPennationModel().calcTendonVelocity(mli.cosPennationAngle,
                mli.sinPennationAngle, dphidt, mli.fiberLength, dlce, dmcldt);
        }

        // Check to see whether the fiber state is clamped.
        double fiberStateClamped = 0.0;
        if(isFiberStateClamped(mli.fiberLength,dlce)) {
            dlce = 0.0;
            dlceN = 0.0;
            dlceAT = 0.0;
            dphidt = 0.0;
            dtl = dmcldt;
            fv = 1.0; //to be consistent with a fiber velocity of 0
            fiberStateClamped = 1.0;
        }

        // Populate the struct.
        fvi.fiberVelocity                = dlce;
        fvi.normFiberVelocity            = dlceN;
        fvi.fiberVelocityAlongTendon     = dlceAT;
        fvi.pennationAngularVelocity     = dphidt;
        fvi.tendonVelocity               = dtl;
        fvi.normTendonVelocity           = dtl/getTendonSlackLength();
        fvi.fiberForceVelocityMultiplier = fv;

        fvi.userDefinedVelocityExtras.resize(1);
        fvi.userDefinedVelocityExtras[0] = fiberStateClamped;

    } catch(const std::exception &x) {
        std::string msg = "Exception caught in Millard2012Muscle::"
                          "calcFiberVelocityInfo from " + getName() + "\n"
                           + x.what();
        throw OpenSim::Exception(msg);
    }
}

//==============================================================================
// MUSCLE INTERFACE REQUIREMENTS -- MUSCLE DYNAMICS INFO
//==============================================================================
void Millard2012Muscle::
calcMuscleDynamicsInfo(const SimTK::State& s, MuscleDynamicsInfo& mdi) const
{
    try {
        // Get the quantities that we've already computed.
        const MuscleLengthInfo &mli = getMuscleLengthInfo(s);
        const FiberVelocityInfo &mvi = getFiberVelocityInfo(s);
        double fiberStateClamped = mvi.userDefinedVelocityExtras[0];

        // Get the properties of this muscle.
        double tendonSlackLen = getTendonSlackLength();
        double optFiberLen    = getOptimalFiberLength();
        double fiso           = getMaxIsometricForce();
        //double penHeight      = penMdl.getParallelogramHeight();
        const TendonForceLengthCurve& fseCurve = get_TendonForceLengthCurve();

        // Compute dynamic quantities.
        double a = SimTK::NaN;
        if(!get_ignore_activation_dynamics()) {
            a = getActivationModel().clampActivation(
                    getActivationModel().calculateActivation(
                            getStateVariableValue(
                                    s, STATE_CALCIUM_CONCENTRATION),
                            mli.fiberLength));
        } else {
            a = getActivationModel().clampActivation(getControl(s));
        }

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
        std::string msg = "Exception caught in Millard2012Muscle::"
                          "calcMuscleDynamicsInfo from " + getName() + "\n"
                          + x.what();
        cerr << msg << endl;
        throw OpenSim::Exception(msg);
    }
}

//==============================================================================
// MODELCOMPONENT INTERFACE REQUIREMENTS
//==============================================================================
void Millard2012Muscle::extendConnectToModel(Model& model)
{
    Super::extendConnectToModel(model);
}

void Millard2012Muscle::
extendAddToSystem(SimTK::MultibodySystem& system) const
{
    Super::extendAddToSystem(system);

    if(!get_ignore_activation_dynamics()) {
        addStateVariable(STATE_CALCIUM_CONCENTRATION);
    }
    if(!get_ignore_tendon_compliance()) {
        addStateVariable(STATE_FIBER_LENGTH_NAME);
    }
}

void Millard2012Muscle::
extendInitStateFromProperties(SimTK::State& s) const
{
    Super::extendInitStateFromProperties(s);

    if(!get_ignore_activation_dynamics()) {
        setActivation(s, getDefaultCalciumConcentration());
    }
    if(!get_ignore_tendon_compliance()) {
        setFiberLength(s, getDefaultFiberLength());
    }
}

void Millard2012Muscle::
extendSetPropertiesFromState(const SimTK::State& s)
{
    Super::extendSetPropertiesFromState(s);

    if(!get_ignore_activation_dynamics()) {
        setDefaultCalciumConcentration(
                getStateVariableValue(s, STATE_CALCIUM_CONCENTRATION));
    }
    if(!get_ignore_tendon_compliance()) {
        setDefaultFiberLength(getStateVariableValue(s,STATE_FIBER_LENGTH_NAME));
    }
}

void Millard2012Muscle::
    computeStateVariableDerivatives(const SimTK::State& s) const
{
    // Activation dynamics if not ignored
    if(!get_ignore_activation_dynamics()) {
        double gammadot = 0;
        // if not disabled or overridden then compute its derivative
        if (appliesForce(s) && !isActuationOverridden(s)) {
            gammadot = getCalciumConcentrationDerivative(s);
        }
        setStateVariableDerivativeValue(
                s, STATE_CALCIUM_CONCENTRATION, gammadot);
    }

    // Fiber length is the next state (if it is a state at all)
    if(!get_ignore_tendon_compliance()) {
        double ldot = 0;
        // if not disabled or overridden then compute its derivative
        if (appliesForce(s) && !isActuationOverridden(s)) {
            ldot = getFiberVelocity(s);
        }
        setStateVariableDerivativeValue(s, STATE_FIBER_LENGTH_NAME, ldot);
    }
}

//==============================================================================
// PRIVATE METHODS
//==============================================================================
SimTK::Vec3 Millard2012Muscle::
calcDampedNormFiberVelocity(double fiso,
                            double a,
                            double fal,
                            double fpe,
                            double fse,
                            double beta,
                            double cosPhi) const
{
    SimTK::Vec4 fiberForceV;
    SimTK::Vec3 result;

    int maxIter = 20; //this routine converges quickly; 20 is quite generous
    double tol = 1.0e-10*fiso;
    if(tol < SimTK::SignificantReal*100) {
        tol = SimTK::SignificantReal*100;
    }
    double perturbation   = 0.0;
    double fiberForce     = 0.0;
    double err            = 1.0e10;
    double derr_d_dlceNdt = 0.0;
    double delta          = 0.0;
    double iter           = 0.0;

    // Get a really excellent starting position to reduce the number of
    // iterations. This reduces the simulation time by about 1%.
    double fv = calcFv(max(a,0.01), max(fal,0.01), fpe, fse, max(cosPhi,0.01));
    double dlceN_dt = fvInvCurve.calcValue(fv);

    // The approximation is poor beyond the maximum velocities.
    if(dlceN_dt > 1.0) {
        dlceN_dt = 1.0;
    }
    if(dlceN_dt < -1.0) {
        dlceN_dt = -1.0;
    }

    double df_d_dlceNdt = 0.0;

    while(abs(err) > tol && iter < maxIter) {
        fv = get_ForceVelocityCurve().calcValue(dlceN_dt);
        fiberForceV = calcFiberForce(fiso,a,fal,fv,fpe,dlceN_dt);
        fiberForce = fiberForceV[0];

        err = fiberForce*cosPhi - fse*fiso;
        df_d_dlceNdt = calc_DFiberForce_DNormFiberVelocity(fiso,a,fal,
                                                           beta,dlceN_dt);
        derr_d_dlceNdt = df_d_dlceNdt*cosPhi;

        if(abs(err) > tol && abs(derr_d_dlceNdt) > SimTK::SignificantReal) {
            delta = -err/derr_d_dlceNdt;
            dlceN_dt = dlceN_dt + delta;

        } else if(abs(derr_d_dlceNdt) < SimTK::SignificantReal) {
            // Perturb the solution if we've lost rank. This should never happen
            // for this problem since dfv_d_dlceNdt > 0 and b > 0 (and so
            // derr_d_dlceNdt > 0).
            perturbation = 2.0*((double)rand())/((double)RAND_MAX)-1.0;
            dlceN_dt = dlceN_dt + perturbation*0.05;
        }
        iter++;
    }

    double converged = 1.0;

    // If we failed to converge, it's probably because the fiber is at its lower
    // bound. That decision is made further down the line, so if convergence
    // didn't happen, let the user know and return a NaN.
    if(abs(err) > tol) {
        dlceN_dt  = -1.0;
        converged = 0.0;
    }

    result[0] = dlceN_dt;
    result[1] = err;
    result[2] = converged;
    return result;
}

double Millard2012Muscle::calcFv(double a,
                                            double fal,
                                            double fp,
                                            double fse,
                                            double cosphi) const
{   return ( fse/cosphi - fp ) / (a*fal); }

SimTK::Vec4 Millard2012Muscle::
calcFiberForce(double fiso,
               double a,
               double fal,
               double fv,
               double fpe,
               double dlceN) const
{
    double beta = getFiberDamping();
    double fa   = fiso * (a*fal*fv);
    double fp1  = fiso * fpe;
    double fp2  = fiso * beta*dlceN;
    double fm   = fa + (fp1+fp2);

    SimTK::Vec4 fiberF;
    fiberF[0] = fm;
    fiberF[1] = fa;
    fiberF[2] = fp1; //conservative passive force
    fiberF[3] = fp2; //non-conservative passive force
    return fiberF;
}

double Millard2012Muscle::calcActivation(double fiso,
                                                    double ftendon,
                                                    double cosPhi,
                                                    double fal,
                                                    double fv,
                                                    double fpe,
                                                    double dlceN) const
{
    double beta = getFiberDamping();
    double activation = 0.0;

    // If the fiber cannot generate any force due to its pennation angle,
    // active-force-length or force-velocity values, leave activation as 0.
    if(cosPhi > SimTK::SignificantReal && fal*fv > SimTK::SignificantReal) {
        activation = ( (ftendon /(fiso*cosPhi)) - fpe - beta*dlceN ) / (fal*fv);
    }
    return activation;
}

double Millard2012Muscle::calcFiberStiffness(double fiso,
                                                        double a,
                                                        double fv,
                                                        double lceN,
                                                        double optFibLen) const
{
    const FiberForceLengthCurve& fpeCurve  = get_FiberForceLengthCurve();
    const ActiveForceLengthCurve& falCurve = get_ActiveForceLengthCurve();
    double DlceN_Dlce = 1.0/optFibLen;
    double Dfal_Dlce  = falCurve.calcDerivative(lceN,1) * DlceN_Dlce;
    double Dfpe_Dlce  = fpeCurve.calcDerivative(lceN,1) * DlceN_Dlce;

    // DFm_Dlce
    return  fiso * (a*Dfal_Dlce*fv + Dfpe_Dlce);
}

double Millard2012Muscle::
calc_DFiberForce_DNormFiberVelocity(double fiso,
                                    double a,
                                    double fal,
                                    double beta,
                                    double dlceN_dt) const
{
    // dfm_d_dlceNdt
    return fiso * (a*fal*get_ForceVelocityCurve().calcDerivative(dlceN_dt,1)
                   + beta);
}

double Millard2012Muscle::
calc_DFiberForceAT_DFiberLength(double fiberForce,
                                double fiberStiffness,
                                double lce,
                                double sinPhi,
                                double cosPhi) const
{
    double Dphi_Dlce    = getPennationModel().
                              calc_DPennationAngle_DfiberLength(lce);
    double Dcosphi_Dlce = -sinPhi*Dphi_Dlce;

    // The stiffness of the fiber along the direction of the tendon. For small
    // changes in length parallel to the fiber, this quantity is
    // D(FiberForceAlongTendon) / D(fiberLength)
    // dFmAT/dlce = d/dlce( fiso * (a *fal*fv + fpe + beta*dlceN)*cosPhi )
    return fiberStiffness*cosPhi + fiberForce*Dcosphi_Dlce;
}

double Millard2012Muscle::
calc_DFiberForceAT_DFiberLengthAT(double dFmAT_d_lce,
                                  double sinPhi,
                                  double cosPhi,
                                  double lce) const
{
    double dphi_d_lce = getPennationModel().
                            calc_DPennationAngle_DfiberLength(lce);

    // The change in length of the fiber length along the tendon.
    // lceAT = lce*cos(phi)
    double DlceAT_Dlce = cosPhi - lce*sinPhi*dphi_d_lce;

    // dFmAT/dlceAT = (dFmAT/dlce)*(1/(dlceAT/dlce))
    //              = dFmAT/dlceAT
    return dFmAT_d_lce * (1.0/DlceAT_Dlce);
}

double Millard2012Muscle::
calc_DTendonForce_DFiberLength(double dFt_d_tl,
                               double lce,
                               double sinphi,
                               double cosphi) const
{
    double dphi_d_lce = getPennationModel().
                            calc_DPennationAngle_DfiberLength(lce);
    double dtl_d_lce  = getPennationModel().
                            calc_DTendonLength_DfiberLength(lce, sinphi, cosphi,
                                                            dphi_d_lce);
    // dFt_d_lce
    return dFt_d_tl*dtl_d_lce;
}

//==============================================================================
// PRIVATE UTILITY CLASS MEMBERS
//==============================================================================
bool Millard2012Muscle::
isFiberStateClamped(double lce, double dlceN) const
{
    bool clamped = false;

    // Get the minimum active fiber length (in meters).
    double minFiberLength = getMinimumFiberLength();

    // If the fiber is clamped and shortening, then the fiber is either shorter
    // than the pennation model allows or shorter than the active-force-length
    // curve allows.
    if( (lce <= minFiberLength && dlceN <= 0) || lce < minFiberLength) {
        clamped = true;
    }
    return clamped;
}

double Millard2012Muscle::clampFiberLength(double lce) const
{   
    return max(lce, getMinimumFiberLength());
}

std::pair<Millard2012Muscle::StatusFromEstimateMuscleFiberState,
          Millard2012Muscle::ValuesFromEstimateMuscleFiberState>
Millard2012Muscle::estimateMuscleFiberState(
                                    const double aActivation,
                                    const double pathLength,
                                    const double pathLengtheningSpeed,
                                    const double aSolTolerance,
                                    const int aMaxIterations,
                                    bool staticSolution) const
{
    // If seeking a static solution, set velocities to zero and avoid the
    // velocity-sharing algorithm below, as it can produce nonzero fiber and
    // tendon velocities even if pathLengtheningSpeed is zero.
    if(abs(pathLengtheningSpeed) < SimTK::SignificantReal) {
        staticSolution = true;
    }

    // Using short variable names to facilitate writing out long equations
    const double ma        = aActivation;
    const double ml        = pathLength;
    const double dml       = pathLengtheningSpeed;
    //Shorter version of the constants
    const double tsl       = getTendonSlackLength();
    const double ofl       = getOptimalFiberLength();
    const double fiso      = getMaxIsometricForce();
    const double vmax      = getMaxContractionVelocity();

    const TendonForceLengthCurve& fseCurve = get_TendonForceLengthCurve();
    const FiberForceLengthCurve& fpeCurve  = get_FiberForceLengthCurve();
    const ActiveForceLengthCurve& falCurve = get_ActiveForceLengthCurve();
    double fse = 0.0;  // normalized tendon (series element) force
    double fal = 0.0;  // normalized active-force-length multiplier
    double fv  = 0.0;  // normalized force-velocity multiplier
    double fpe = 0.0;  // normalized parallel element force

    // Position level
    double tl  = getTendonSlackLength()*1.01;  // begin with small tendon force
    double lce = clampFiberLength(getPennationModel().calcFiberLength(ml,tl));

    double phi = 0.0;
    double cosphi = 1.0;
    double sinphi = 0.0;

    //Normalized quantities
    double tlN    = tl/tsl;
    double lceN   = lce/ofl;
    // Velocity level first guess assume static
    double dtl    = 0.0;
    double dlce = 0.0;
    double dlceN  = 0.0;

    // Internal variables for the loop
    double Fm           = 0.0;   // fiber force
    double FmAT         = 0.0;   // fiber force along tendon
    double Ft           = 0.0;   // tendon force
    double ferr         = SimTK::MostPositiveReal;  // solution error

    double dFm_dlce     = 0.0;   // partial of muscle force w.r.t. lce
    double dFmAT_dlce   = 0.0;   // partial of muscle force along tl w.r.t. lce
    double dFmAT_dlceAT = 0.0;   // partial of muscle force along tl w.r.t. lce
                                 //     along the tendon
    double dFt_d_lce    = 0.0;   // partial of tendon force w.r.t. lce
    double dFt_d_tl     = 0.0;   // partial of tendon force w.r.t. tl
    double dferr_d_lce  = 0.0;   // partial of solution error w.r.t lce
    double delta_lce    = 0.0;   // change in lce

    SimTK::Vec4 fiberForceV(SimTK::NaN);

    //*******************************
    // Helper functions
    //Update position level quantities, only if they won't go singular
    auto positionFunc = [&] {
        phi = getPennationModel().calcPennationAngle(lce);
        cosphi = cos(phi);
        sinphi = sin(phi);
        tl = ml - lce*cosphi;
        lceN = lce / ofl;
        tlN = tl / tsl;
    };

    // Functional to update the force multipliers
    auto multipliersFunc = [&] {
        fal = falCurve.calcValue(lceN);
        fpe = fpeCurve.calcValue(lceN);
        fse = fseCurve.calcValue(tlN);
    };

    // Functional to compute the equilibrium force error
    auto ferrFunc = [&] {
        fiberForceV = calcFiberForce(fiso, ma, fal, fv, fpe, dlceN);
        Fm = fiberForceV[0];
        FmAT = Fm * cosphi;
        Ft = fse*fiso;
        ferr = FmAT - Ft;
    };

    // Functional to compute the partial derivative of muscle force w.r.t. lce
    auto partialsFunc = [&] {
        dFm_dlce = calcFiberStiffness(fiso, ma, fv, lceN, ofl);
        dFmAT_dlce = calc_DFiberForceAT_DFiberLength(Fm, dFm_dlce, lce,
            sinphi, cosphi);
        dFmAT_dlceAT = calc_DFiberForceAT_DFiberLengthAT(dFmAT_dlce, sinphi,
            cosphi, lce);
        dFt_d_tl = fseCurve.calcDerivative(tlN, 1)*fiso / tsl;
        dFt_d_lce = calc_DTendonForce_DFiberLength(dFt_d_tl, lce,
            sinphi, cosphi);
    };

    // Functional to estimate fiber velocity and force-velocity multiplier
    // from the relative fiber and tendon stiffnesses from above partials 
    auto velocityFunc = [&] {
        /* Update velocity-level quantities if not staticSolution.
        Share the muscle velocity between the tendon and the fiber
        according to their relative stiffnesses:

        Fm-Ft = 0                     Equilibrium equation   [1]
        d/dt Fm - d/dt Ft = 0         Time derivative        [2]
        lp = lm + lt                  Path definition        [3]
        d/dt lp = d/dt lm + d/dt lt   Path derivative        [4]

        Computing a linearized model of [2]:
        Fm = Fm0 + Km*lceAT                                  [5]
        Ft = Ft0 Kt*xt                                       [6]

        Taking its time derivative:
        dFm_d_xm = Km*dlceAT + dKm_d_t*lceAT (assume dKm_d_t = 0)   [7]
        dFt_d_xt = Kt*dtl + dKt_d_t*dtl (assume dKt_d_t = 0)        [8]

        Substituting 7 and 8 into 2:
        Km dlceAT - Kt dtl = 0

        Using Eqn 4, we have 2 equations in 2 unknowns. Can now solve for
        tendon velocity, or the velocity of the fiber along the tendon.

        This is a heuristic. The above assumptions are necessary since
        computing the partial derivatives of Km or Kt requires acceleration-
        level knowledge, which is not available in general.

        Stiffness of the muscle is the stiffness of the tendon and the fiber
        (along the tendon) in series. */

        if (!staticSolution) {
            // The "if" statement here is to handle the special case where the
            // negative stiffness of the fiber (which happens in this model) is
            // equal to the positive stiffness of the tendon.
            if (abs(dFmAT_dlceAT + dFt_d_tl) > SimTK::SignificantReal
                && tlN > 1.0) {
                dtl = dFmAT_dlceAT / (dFmAT_dlceAT + dFt_d_tl) * dml;
            }
            else {
                dtl = dml;
            }

            // Update fiber velocity
            dlce = getPennationModel().calcFiberVelocity(cosphi, dml, dtl);
            dlceN = dlce / (vmax*ofl);
            // Update the force-velocity multiplier
            fv = get_ForceVelocityCurve().calcValue(dlceN);
        }
    };

    //*******************************
    //Initialize the loop
    int iter = 0;

    // Estimate the position level quantities (lengths, angles) of the muscle
    positionFunc();

    // Multipliers based on initial fiber-length estimate
    multipliersFunc();

    // Starting guess at the force-velocity multiplier is static
    fv = 1.0;

    fiberForceV = calcFiberForce(fiso, ma, fal, fv, fpe, dlceN);
    Fm = fiberForceV[0];

    // Compute the partial derivative of the force error w.r.t. lce
    partialsFunc();

    // update fiber and tendon velocity and velocity multiplier
    velocityFunc();

    // Compute the force error
    ferrFunc();

    // Update the partial derivatives of the force error w.r.t. lce with 
    // newly estimated fv
    partialsFunc();

    double ferrPrev = ferr;
    double lcePrev = lce;

    double h =1.0;
    while( (abs(ferr) > aSolTolerance) && (iter < aMaxIterations)) {
        // Compute the search direction
        dferr_d_lce = dFmAT_dlce - dFt_d_lce;
        h = 1.0;

        while (abs(ferr) >= abs(ferrPrev)) {
            // Compute the Newton step
            delta_lce = -h*ferrPrev / dferr_d_lce;
            // Take a Newton Step if the step is nonzero
            if (abs(delta_lce) > SimTK::SignificantReal)
                lce = lcePrev + delta_lce;
            else {
                // We've stagnated or hit a limit; assume we are hitting local
                // minimum and attempt to approach from the other direction.
                lce = lcePrev - sign(delta_lce)*SimTK::SqrtEps;
                // Force a break, which will update the derivatives of
                // the muscle force and estimate of the fiber-velocity
                h = 0;
            }

            if (lce < getMinimumFiberLength()) {
                lce = getMinimumFiberLength();
            }

            // Update the muscles's position level quantities (lengths, angles)
            positionFunc();

            // Update the muscle force multipliers
            multipliersFunc();

            // Compute the force error assuming fiber-velocity is unchanged
            ferrFunc();

            if (h <= SimTK::SqrtEps) {
                break;
            }
            else
                h = 0.5*h;
        }

        ferrPrev = ferr;
        lcePrev = lce;
        
        // Update the partial derivative of the force error w.r.t. lce
        partialsFunc();
        // Update velocity estimate and velocity multiplier
        velocityFunc();

        iter++;
    }

    // Populate the result map.
    ValuesFromEstimateMuscleFiberState resultValues;

    if(abs(ferr) < aSolTolerance) {  // The solution converged.

        if (isFiberStateClamped(lce, dlceN)) {
            lce = getMinimumFiberLength();
        }

        resultValues["solution_error"] = ferr;
        resultValues["iterations"]     = (double)iter;
        resultValues["fiber_length"]   = lce;
        resultValues["fiber_velocity"] = dlce;
        resultValues["tendon_force"]   = fse*fiso;

        return std::pair<StatusFromEstimateMuscleFiberState,
                         ValuesFromEstimateMuscleFiberState>
          (StatusFromEstimateMuscleFiberState::Success_Converged, resultValues);
    }

    // Fiber length is at or exceeds its lower bound.
    if (lce <= getMinimumFiberLength()) {

        lce = getMinimumFiberLength();
        phi    = getPennationModel().calcPennationAngle(lce);
        cosphi = cos(phi);
        tl     = getPennationModel().calcTendonLength(cosphi,lce,ml);
        lceN   = lce/ofl;
        tlN    = tl/tsl;
        fse    = fseCurve.calcValue(tlN);

        resultValues["solution_error"] = ferr;
        resultValues["iterations"]     = (double)iter;
        resultValues["fiber_length"]   = lce;
        resultValues["fiber_velocity"] = 0;
        resultValues["tendon_force"]   = fse*fiso;

        return std::pair<StatusFromEstimateMuscleFiberState,
                         ValuesFromEstimateMuscleFiberState>
            (StatusFromEstimateMuscleFiberState::Warning_FiberAtLowerBound,
             resultValues);
    }

    resultValues["solution_error"] = ferr;
    resultValues["iterations"]     = (double)iter;
    resultValues["fiber_length"]   = SimTK::NaN;
    resultValues["fiber_velocity"] = SimTK::NaN;
    resultValues["tendon_force"]   = SimTK::NaN;

    return std::pair<StatusFromEstimateMuscleFiberState,
                        ValuesFromEstimateMuscleFiberState>
        (StatusFromEstimateMuscleFiberState::Failure_MaxIterationsReached,
            resultValues);
}