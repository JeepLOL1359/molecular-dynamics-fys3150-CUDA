#include "system.h"
#include "velocityverlet.h"
#include "lennardjones.h"
#include "statisticssampler.h"
#include "unitconverter.h"
#include "math/random.h"

System::System()
{

}

System::~System()
{
    for(Atom *atom : m_atoms) {
        delete atom;
    }
    m_atoms.clear();
}

void System::applyPeriodicBoundaryConditions() {
    // Read here: http://en.wikipedia.org/wiki/Periodic_boundary_conditions#Practical_implementation:_continuity_and_the_minimum_image_convention
}

void System::removeTotalMomentum() {
    // Find the total momentum and remove momentum equally on each atom so the total momentum becomes zero.
}

void System::createFCCLattice(int numberOfUnitCellsEachDimension,
    double latticeConstant,
    double temperature)
{
    m_atoms.clear();

    int n = numberOfUnitCellsEachDimension;

    double mass = UnitConverter::massFromSI(6.63352088e-26);

    // FCC offsets inside one unit cell
    vec3 offsets[4] = {
        vec3(0.0, 0.0, 0.0),
        vec3(0.5, 0.5, 0.0),
        vec3(0.5, 0.0, 0.5),
        vec3(0.0, 0.5, 0.5)
    };

    for (int x = 0; x < n; x++) {
        for (int y = 0; y < n; y++) {
            for (int z = 0; z < n; z++) {

                vec3 base = vec3(x, y, z) * latticeConstant;

                for (int i = 0; i < 4; i++) {

                    Atom* atom = new Atom(mass);

                    atom->position = base + offsets[i] * latticeConstant;

                    atom->resetVelocityMaxwellian(temperature);

                    m_atoms.push_back(atom);
                }
            }
        }
    }

    // Box size = n * latticeConstant
    setSystemSize(vec3(n, n, n) * latticeConstant);
}

void System::calculateForces() {
    for(Atom *atom : m_atoms) {
        atom->resetForce();
    }
    m_potential.calculateForces(*this); // this is a pointer, *this is a reference to this object
}

void System::step(double dt) {
    m_integrator.integrate(*this, dt);
    m_steps++;
    m_time += dt;
}
