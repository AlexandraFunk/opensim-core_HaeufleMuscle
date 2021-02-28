#include "OpenSim/Common/STOFileAdapter.h"
#include "OpenSim/Common/TimeSeriesTable.h"

#include <OpenSim/Common/IO.h>
#include <OpenSim/OpenSim.h> 
using namespace SimTK;
using namespace OpenSim;

int main() {
    Model model;
    model.setName("bicep_curl");
    model.setUseVisualizer(true);

    // Create two links, each with a mass of 1 kg, center of mass at the body's
    // origin, and moments and products of inertia of zero.
    OpenSim::Body* humerus = new OpenSim::Body("humerus", 1, Vec3(0), Inertia(0));
    OpenSim::Body* radius  = new OpenSim::Body("radius",  1, Vec3(0), Inertia(0));

    // Connect the bodies with pin joints. Assume each body is 1 m long.
    PinJoint* shoulder = new PinJoint("shoulder",
            // Parent body, location in parent, orientation in parent.
            model.getGround(), Vec3(0), Vec3(0),
            // Child body, location in child, orientation in child.
            *humerus, Vec3(0, 1, 0), Vec3(0));
    PinJoint* elbow = new PinJoint("elbow",
            *humerus, Vec3(0), Vec3(0), *radius, Vec3(0, 1, 0), Vec3(0));

    // Add a muscle that flexes the elbow.
    // Millard2012EquilibriumMuscle* biceps = new
    //    Millard2012EquilibriumMuscle("biceps", 200, 0.6, 0.55, 0);

    //Haeufle2014Muscle* biceps =
    //        new Haeufle2014Muscle("biceps", 200, 0.45, 0.25, 0.0);
    //Haeufle2014Muscle* biceps = new Haeufle2014Muscle(
    //        "biceps", 200, 0.6, 0.55, 20.0 * SimTK_DEGREE_TO_RADIAN);
 
    Haeufle2014Muscle* biceps = new Haeufle2014Muscle(
            "biceps", 200, 0.6, 0.15, 10.0 * SimTK_DEGREE_TO_RADIAN);
    biceps->set_fibre_maximum_damping_dpe(0.3);
    biceps->set_fibre_offset_damping_rpe(0.01);

    // set damping params to 0 to be more similar to millard muscle
    // biceps->setParallelDampingParams(0.0, 0.0);
    // biceps->setTendonDampingParams(0.0, 0.0);
    // printf("Dse: %f\n", biceps->get_dse_damping_factor());
    // printf("Rse: %f\n", biceps->get_rse_damping_factor());
    // printf("Dpe: %f\n", biceps->get_dpe_damping_factor());
    // printf("Rpe: %f\n", biceps->get_rpe_damping_factor());
    // biceps->set_default_calcium_concentration(0.05);


    biceps->addNewPathPoint("origin",    *humerus, Vec3(0, 0.8, 0));
    biceps->addNewPathPoint("insertion", *radius,  Vec3(0, 0.7, 0));

    // Add a controller that specifies the excitation of the muscle.
    PrescribedController* brain = new PrescribedController();
    brain->addActuator(*biceps);
    // Muscle excitation is 0.3 for the first 0.5 seconds, then increases to 1.
    double x[] = {0.0, 0.5, 1.0, 2.0, 5.0};
    double y[] = {0.3, 0.3, 1.0, 0.2, 0.2};
    auto control_function = new PiecewiseLinearFunction(5, x, y);
    brain->prescribeControlForActuator("biceps", control_function);
    //brain->prescribeControlForActuator(
    //        "biceps", new StepFunction(0.5, 3, 0.3, 1.0));

    // Add components to the model.
    model.addBody(humerus);    model.addBody(radius);
    model.addJoint(shoulder);  model.addJoint(elbow);
    model.addForce(biceps);
    model.addController(brain);

    // Add a console reporter to print the muscle fiber force and elbow angle.
    ConsoleReporter* reporter = new ConsoleReporter();
    TableReporter* length_reporter = new TableReporter();
    length_reporter->set_report_time_interval(0.001);
    length_reporter->addToReport(biceps->getOutput("fiber_length"));
    length_reporter->addToReport(biceps->getOutput("tendon_length"));
    length_reporter->addToReport(biceps->getOutput("pennation_angle"));
    length_reporter->addToReport(biceps->getOutput("activation"));
    length_reporter->addToReport(biceps->getOutput("excitation"));
    length_reporter->addToReport(biceps->getOutput("fiber_velocity"));
    length_reporter->addToReport(biceps->getOutput("speed"));


    TableReporter* force_reporter = new TableReporter();
    force_reporter->setName("Forces");
    force_reporter->set_report_time_interval(0.001);
    force_reporter->addToReport(biceps->getOutput("tendon_force"));
    force_reporter->addToReport(biceps->getOutput("active_fiber_force"));
    force_reporter->addToReport(biceps->getOutput("Fpee"));
    force_reporter->addToReport(biceps->getOutput("Fpde"));
    force_reporter->addToReport(biceps->getOutput("Fsee"));
    force_reporter->addToReport(biceps->getOutput("Fsde"));

    model.addComponent(length_reporter);
    model.addComponent(force_reporter);


    // Add display geometry.
    Ellipsoid bodyGeometry(0.1, 0.5, 0.1);
    bodyGeometry.setColor(Gray);
    // Attach an ellipsoid to a frame located at the center of each body.
    PhysicalOffsetFrame* humerusCenter = new PhysicalOffsetFrame(
        "humerusCenter", *humerus, Transform(Vec3(0, 0.5, 0)));
    humerus->addComponent(humerusCenter);
    humerusCenter->attachGeometry(bodyGeometry.clone());
    PhysicalOffsetFrame* radiusCenter = new PhysicalOffsetFrame(
        "radiusCenter", *radius, Transform(Vec3(0, 0.5, 0)));
    radius->addComponent(radiusCenter);
    radiusCenter->attachGeometry(bodyGeometry.clone());

    // Configure the model.
    State& state = model.initSystem();

    Manager manager(model);
    manager.setIntegratorMethod(Manager::IntegratorMethod::RungeKuttaMerson);
    manager.setIntegratorAccuracy(10e-6);

    // Fix the shoulder at its default angle and begin with the elbow flexed.
    shoulder->getCoordinate().setLocked(state, true);
    elbow->getCoordinate().setValue(state, 0.5 * Pi);
    //elbow->getCoordinate().setValue(state, 0.001 * Pi);
    model.equilibrateMuscles(state);

    // Configure the visualizer.
    model.updMatterSubsystem().setShowDefaultGeometry(true);
    Visualizer& viz = model.updVisualizer().updSimbodyVisualizer();
    viz.setBackgroundType(viz.SolidColor);
    viz.setBackgroundColor(White);

    // Simulate.
    //simulate(model, state, 5.0);

    state.setTime(0.0);
    manager.initialize(state);
    std::cout << "\nIntegrating from " << 0.0 << " to " << 5.0
              << std::endl;
    manager.integrate(5.0);


    // Display the table reporters
    auto forcesTable = force_reporter->getTable();
    STOFileAdapter_<double>::write(forcesTable, "forces_states.sto");
    auto lengthTable = length_reporter->getTable();
    STOFileAdapter_<double>::write(lengthTable, "length_states.sto");

    // To print (serialize) the latest connections of the model, it is
    // necessary to finalizeConnections() first.
    model.finalizeConnections();
    // Save the OpenSim model to a file
    model.print("bicepsMain.osim");

    std::cout << "\nFinished successfully!" << std::endl;

    return 0;
};