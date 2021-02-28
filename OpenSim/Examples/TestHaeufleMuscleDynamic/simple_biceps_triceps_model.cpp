#include <OpenSim/OpenSim.h> 
using namespace SimTK;
using namespace OpenSim;

int main() {
    Model model;
    model.setName("bicep_trizeps_curl");
    model.setUseVisualizer(true);

    // Create two links, each with a mass of 1 kg, center of mass at the body's
    // origin, and moments and products of inertia of zero.
    OpenSim::Body* humerus =
            new OpenSim::Body("humerus", 1, Vec3(0), Inertia(0));
    OpenSim::Body* radius = new OpenSim::Body("radius", 1, Vec3(0), Inertia(0));

    // Connect the bodies with pin joints. Assume each body is 1 m long.
    PinJoint* shoulder = new PinJoint("shoulder",
            // Parent body, location in parent, orientation in parent.
            model.getGround(), Vec3(0), Vec3(0),
            // Child body, location in child, orientation in child.
            *humerus, Vec3(0, 1, 0), Vec3(0));
    PinJoint* elbow = new PinJoint("elbow", *humerus, Vec3(0), Vec3(0), *radius,
            Vec3(0, 1, 0), Vec3(0));

    // Add a muscle that flexes the elbow.
    // Millard2012EquilibriumMuscle* biceps = new
    //    Millard2012EquilibriumMuscle("biceps", 200, 0.6, 0.55, 0);

    Haeufle2014Muscle* biceps =
            new Haeufle2014Muscle("biceps", 200, 0.45, 0.25, 0.0);

    Haeufle2014Muscle* triceps =
            new Haeufle2014Muscle("triceps", 200, 0.9, 0.5, 0.0);

    // set damping params to 0 to be more similar to millard muscle
    // biceps->setParallelDampingParams(0.0, 0.0);
    // biceps->setTendonDampingParams(0.0, 0.0);
    // printf("Dse: %f\n", biceps->get_dse_damping_factor());
    // printf("Rse: %f\n", biceps->get_rse_damping_factor());
    // printf("Dpe: %f\n", biceps->get_dpe_damping_factor());
    // printf("Rpe: %f\n", biceps->get_rpe_damping_factor());
    // biceps->set_default_calcium_concentration(0.05);

    biceps->addNewPathPoint("origin", *humerus, Vec3(0, 0.8, 0));
    biceps->addNewPathPoint("insertion", *radius, Vec3(0, 0.7, 0));
    triceps->addNewPathPoint("origin", *humerus, Vec3(-0.2, 0.8, 0));
    triceps->addNewPathPoint("ViaPoint1", *humerus, Vec3(0, -0.3, 0));
    triceps->addNewPathPoint("insertion", *radius, Vec3(0, 0.7, 0));

    // Add a controller that specifies the excitation of the muscle.
    PrescribedController* brain = new PrescribedController();
    brain->addActuator(*biceps);
    brain->addActuator(*triceps);

    // creating hat function of type
    // t = [0,10]
    // x = 0 for 1 sec then 1 for 1 sec
    double x[10] = {1, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    double y[10] = {1, 0, 1, 0, 1, 0, 1, 0, 1, 0};
    double y_flipped[10] = {0, 1, 0, 1, 0, 1, 0, 1, 0, 1};
    auto controlFunction = new PiecewiseConstantFunction(10, x, y);
    auto controlFunction_flipped =
            new PiecewiseConstantFunction(10, x, y_flipped);
    
    //new StepFunction(0, 1, 0.1, 1)

    brain->prescribeControlForActuator("biceps", controlFunction);
    brain->prescribeControlForActuator("triceps", controlFunction_flipped);
    
    // Add components to the model.
    model.addBody(humerus);    model.addBody(radius);
    model.addJoint(shoulder);  model.addJoint(elbow);
    model.addForce(biceps);
    model.addForce(triceps);
    model.addController(brain);

    /*
    // Add a console reporter to print the muscle fiber force and elbow angle.
    ConsoleReporter* reporter = new ConsoleReporter();
    reporter->set_report_time_interval(1.0);
    reporter->addToReport(biceps->getOutput("fiber_force"));
    reporter->addToReport(biceps->getOutput("activation"));
    reporter->addToReport(
        elbow->getCoordinate(PinJoint::Coord::RotationZ).getOutput("value"),
        "elbow_angle");
    */

    ConsoleReporter* reportTable = new ConsoleReporter();
    reportTable->set_report_time_interval(1);
    reportTable->addToReport(biceps->getOutput("fiber_force"));
    reportTable->addToReport(biceps->getOutput("tendon_force"));
    reportTable->addToReport(biceps->getOutput("activation"));
    reportTable->addToReport(biceps->getOutput("fiber_length"));
    reportTable->addToReport(biceps->getOutput("tendon_length"));
    reportTable->addToReport(triceps->getOutput("fiber_force"));
    reportTable->addToReport(triceps->getOutput("tendon_force"));
    reportTable->addToReport(triceps->getOutput("activation"));
    reportTable->addToReport(triceps->getOutput("fiber_length"));
    reportTable->addToReport(triceps->getOutput("tendon_length"));
    reportTable->addToReport(
            elbow->getCoordinate(PinJoint::Coord::RotationZ).getOutput("value"),
            "elbow_angle");
    reportTable->addToReport(
        shoulder->getCoordinate(PinJoint::Coord::RotationZ).getOutput("value"),
            "shoulder_angle");

    model.addComponent(reportTable);

    // Add display geometry.
    /* 
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
    */
    // Configure the model.
    State& state = model.initSystem();
    // Fix the shoulder at its default angle and begin with the elbow flexed.
    shoulder->getCoordinate().setLocked(state, true);
    elbow->getCoordinate().setValue(state, 0.5 * Pi);
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