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

    SimTK::Vec3 translation1(0.1, 0.1, 0.15);
    SimTK::Vec3 translation2(0.1, 0.8, 0.05);

    SimTK::Vec3 xyz_body_rotation1(45.0 * SimTK::Pi / 180.0,
            0.0, 90.0 * SimTK::Pi / 180.0);
    SimTK::Vec3 xyz_body_rotation2(50.1489 * SimTK::Pi / 180.0,
            1.72794 * SimTK::Pi / 180.0, 80.1489 * SimTK::Pi / 180.0);

    OpenSim::WrapViaEllipse* firstEllipse = new OpenSim::WrapViaEllipse();
    firstEllipse->setName("ViaEllipse01");
    firstEllipse->setSemiAxisLengthG(.05);
    firstEllipse->setSemiAxisLengthH(.025);
    firstEllipse->set_translation(translation1);
    firstEllipse->set_xyz_body_rotation(xyz_body_rotation1);

    
    OpenSim::WrapViaEllipse* secondEllipse = new OpenSim::WrapViaEllipse();
    secondEllipse->setName("ViaEllipse02");
    secondEllipse->setSemiAxisLengthG(.01);
    secondEllipse->setSemiAxisLengthH(.02);
    secondEllipse->set_translation(translation2);
    secondEllipse->set_xyz_body_rotation(xyz_body_rotation2);
    

    humerus->addWrapObject(firstEllipse);
    radius->addWrapObject(secondEllipse);

    // Connect the bodies with pin joints. Assume each body is 1 m long.
    PinJoint* shoulder = new PinJoint("shoulder",
            // Parent body, location in parent, orientation in parent.
            model.getGround(), Vec3(0), Vec3(0),
            // Child body, location in child, orientation in child.
            *humerus, Vec3(0, 1, 0), Vec3(0));
    PinJoint* elbow = new PinJoint("elbow",
            *humerus, Vec3(0), Vec3(0), *radius, Vec3(0, 1, 0), Vec3(0));

    // Add a muscle that flexes the elbow.
    Millard2012EquilibriumMuscle* biceps = new
        Millard2012EquilibriumMuscle("biceps", 200, 0.6, 0.55, 0);
    biceps->addNewPathPoint("origin",    *humerus, Vec3(0, 0.8, 0));
    biceps->addNewPathPoint("insertion", *radius,  Vec3(0, 0.7, 0));

    biceps->updGeometryPath().addPathWrap(*firstEllipse);
    biceps->updGeometryPath().addPathWrap(*secondEllipse);
    

    // Add a controller that specifies the excitation of the muscle.
    PrescribedController* brain = new PrescribedController();
    brain->addActuator(*biceps);
    // Muscle excitation is 0.3 for the first 0.5 seconds, then increases to 1.
    brain->prescribeControlForActuator("biceps",
            new StepFunction(0.5, 3, 0.3, 1));

    // Add components to the model.
    model.addBody(humerus);    model.addBody(radius);
    model.addJoint(shoulder);  model.addJoint(elbow);
    model.addForce(biceps);
    model.addController(brain);

    // Add a console reporter to print the muscle fiber force and elbow angle.
    ConsoleReporter* reporter = new ConsoleReporter();
    reporter->set_report_time_interval(0.001);
    //reporter->addToReport(biceps->getOutput("fiber_force"));
    reporter->addToReport(
        elbow->getCoordinate(PinJoint::Coord::RotationZ).getOutput("value"),
        "elbow_angle");
    reporter->addToReport(biceps->getOutput("tendon_length"));
    reporter->addToReport(biceps->getOutput("fiber_length"));
    reporter->addToReport(firstEllipse->getOutput("angleOnEllipse"));
    reporter->addToReport(secondEllipse->getOutput("angleOnEllipse"));

    /*
    std::cout << "Ellipse Outputs:" << std::endl;
    for (int i = 0; i < firstEllipse->getOutputNames().size(); i++) {
        std::cout << "\t" << firstEllipse->getOutputNames()[i]
                  << std::endl;
    }
    */


    //reporter->addToReport(biceps->getGeometryPath().getWrapSet().get(0).getOutput("angleOnEllipse"));
    model.addComponent(reporter);

    // Add display geometry.
    //Ellipsoid bodyGeometry(0.05, 0.05, 0.0);
    //bodyGeometry.setColor(Gray);
    // Attach an ellipsoid to a frame located at the center of each body.
    PhysicalOffsetFrame* humerusCenter = new PhysicalOffsetFrame(
        "humerusCenter", *humerus, Transform(Vec3(0, 0.5, 0)));
    humerus->addComponent(humerusCenter);
    //humerusCenter->attachGeometry(bodyGeometry.clone());
    PhysicalOffsetFrame* radiusCenter = new PhysicalOffsetFrame(
        "radiusCenter", *radius, Transform(Vec3(0, 0.5, 0)));
    radius->addComponent(radiusCenter);
    //radiusCenter->attachGeometry(bodyGeometry.clone());

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