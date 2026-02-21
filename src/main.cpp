#include "math/random.h"
#include "lennardjones.h"
#include "velocityverlet.h"
#include "system.h"
#include "statisticssampler.h"
#include "atom.h"
#include "io.h"
#include "unitconverter.h"
#include <iostream>
#include <iomanip>

using namespace std;

int main(int numberOfArguments, char **argumentList)
{
    int numberOfUnitCells = 4;
    double initialTemperature = UnitConverter::temperatureFromSI(300.0); // measured in Kelvin
    double latticeConstant = UnitConverter::lengthFromAngstroms(5.26); // measured in angstroms

    // If a first argument is provided, it is the number of unit cells
    if(numberOfArguments > 1) numberOfUnitCells = atoi(argumentList[1]);
    // If a second argument is provided, it is the initial temperature (measured in kelvin)
    if(numberOfArguments > 2) initialTemperature = UnitConverter::temperatureFromSI(atof(argumentList[2]));
    // If a third argument is provided, it is the lattice constant determining the density (measured in angstroms)
    if(numberOfArguments > 3) latticeConstant = UnitConverter::lengthFromAngstroms(atof(argumentList[3]));

    double dt = UnitConverter::timeFromSI(1e-15); // Measured in seconds.

    cout << "========== UNIT SYSTEM ==========\n";
    cout << "Length unit      : " << UnitConverter::lengthToSI(1.0) << " m\n";
    cout << "Velocity unit    : " << UnitConverter::velocityToSI(1.0) << " m/s\n";
    cout << "Time unit        : " << UnitConverter::timeToSI(1.0) << " s\n";
    cout << "Mass unit        : " << UnitConverter::massToSI(1.0) << " kg\n";
    cout << "Temperature unit : " << UnitConverter::temperatureToSI(1.0) << " K\n\n";

    System system;
    system.createFCCLattice(numberOfUnitCells, latticeConstant, initialTemperature);
    system.potential().setEpsilon(1.0);
    system.potential().setSigma(1.0);

    system.removeTotalMomentum();

    cout << "========== SIMULATION ==========\n";
    cout << "Number of atoms  : " << system.atoms().size() << "\n";
    cout << "Time step (dt)   : " << dt << "\n";
    cout << "Total steps      : 1000\n\n";

    StatisticsSampler statisticsSampler;
    IO movie("movie.xyz"); // To write the state to file

    cout << "========== ENERGY EVOLUTION ==========\n";
    cout << setw(20) << "Timestep" <<
            setw(20) << "Time" <<
            setw(20) << "Temperature" <<
            setw(20) << "KineticEnergy" <<
            setw(20) << "PotentialEnergy" <<
            setw(20) << "TotalEnergy" << endl;
    for(int timestep=0; timestep<1000; timestep++) {
        system.step(dt);
        statisticsSampler.sample(system);
        if( timestep % 100 == 0 ) {
            // Print the timestep every 100 timesteps
            cout << setw(20) << system.steps() <<
                    setw(20) << system.time() <<
                    setw(20) << statisticsSampler.temperature() <<
                    setw(20) << statisticsSampler.kineticEnergy() <<
                    setw(20) << statisticsSampler.potentialEnergy() <<
                    setw(20) << statisticsSampler.totalEnergy() << endl;
        }
        movie.saveState(system);
    }

    movie.close();

    return 0;
}
