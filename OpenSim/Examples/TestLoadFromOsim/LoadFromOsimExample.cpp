#include <OpenSim/OpenSim.h>

using namespace SimTK;
using namespace OpenSim;

int main() {
    Model model("dpendulum.osim");
    model.setUseVisualizer(true);

    // get the muscle m_flexor from osim file
    Millard2012EquilibriumMuscle& biceps =
            static_cast<Millard2012EquilibriumMuscle&> (model.getMuscles()[0]);

    biceps.set_max_isometric_force(200.0);
    biceps.set_optimal_fiber_length(0.6);
    biceps.set_tendon_slack_length(0.55);
    biceps.set_pennation_angle_at_optimal(0);

    OpenSim::Body& humerus = model.getBodySet()[0];
    OpenSim::Body& radius = model.getBodySet()[1];

    OpenSim::WrapViaEllipse& firstEllipse =
            static_cast<WrapViaEllipse&> (humerus.getWrapObjectSet()[0]);
    OpenSim::WrapViaEllipse& secondEllipse =
            static_cast<WrapViaEllipse&>(radius.getWrapObjectSet()[0]);



    CustomJoint& shoulder = static_cast<CustomJoint&>(model.getJointSet()[0]);
    CustomJoint& elbow = static_cast<CustomJoint&>(model.getJointSet()[1]);

    // Add a controller that specifies the excitation of the muscle.
    PrescribedController* brain = new PrescribedController();
    brain->addActuator(biceps);
    // Muscle excitation is 0.3 for the first 0.5 seconds, then increases to 1.
    brain->prescribeControlForActuator(
            "M_flexor", new StepFunction(0.5, 3, 0.3, 1));

    model.addController(brain);

        // Add a console reporter to print the muscle fiber force and elbow angle.
    ConsoleReporter* reporter = new ConsoleReporter();
    reporter->set_report_time_interval(0.001);
    reporter->addToReport(
            elbow.getCoordinate(0).getOutput("value"),
            "elbow_angle");
    reporter->addToReport(biceps.getOutput("tendon_length"));
    reporter->addToReport(biceps.getOutput("fiber_length"));
    reporter->addToReport(firstEllipse.getOutput("angleOnEllipse"));
    reporter->addToReport(secondEllipse.getOutput("angleOnEllipse"));


    // reporter->addToReport(biceps->getGeometryPath().getWrapSet().get(0).getOutput("angleOnEllipse"));
    model.addComponent(reporter);



    // Configure the model.
    State& state = model.initSystem();
    // Fix the shoulder at its default angle and begin with the elbow flexed.
    shoulder.getCoordinate().setLocked(state, true);
    elbow.getCoordinate().setValue(state, 0.5 * Pi);
    model.equilibrateMuscles(state);

    // Configure the visualizer.
    model.updMatterSubsystem().setShowDefaultGeometry(true);
    Visualizer& viz = model.updVisualizer().updSimbodyVisualizer();
    viz.setBackgroundType(viz.SolidColor);
    viz.setBackgroundColor(White);

    // Simulate.
    simulate(model, state, 10.0);

    return 0;
};