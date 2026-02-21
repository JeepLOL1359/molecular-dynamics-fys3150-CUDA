#include "lennardjones.h"
#include "system.h"

double LennardJones::potentialEnergy() const
{
    return m_potentialEnergy;
}

double LennardJones::sigma() const
{
    return m_sigma;
}

void LennardJones::setSigma(double sigma)
{
    m_sigma = sigma;
}

double LennardJones::epsilon() const
{
    return m_epsilon;
}

void LennardJones::setEpsilon(double epsilon)
{
    m_epsilon = epsilon;
}

void LennardJones::calculateForces(System& system)
{
    m_potentialEnergy = 0; // Remember to compute this in the loop

    auto& atoms = system.atoms();
    int N = atoms.size();

    double cutoff = 3.0 * m_sigma;
    double cutoff2 = cutoff * cutoff;

    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {

            vec3 rij = atoms[i]->position - atoms[j]->position;

            double r2 = rij.lengthSquared();

            if (r2 > cutoff2) continue;

            double inv_r2 = 1.0 / r2;

            double sigma2 = m_sigma * m_sigma;
            double sigma6 = sigma2 * sigma2 * sigma2;
            double sigma12 = sigma6 * sigma6;

            double inv_r6 = inv_r2 * inv_r2 * inv_r2;
            double inv_r12 = inv_r6 * inv_r6;

            // Potential energy
            double potential = 4.0 * m_epsilon * (sigma12 * inv_r12 - sigma6 * inv_r6);
            m_potentialEnergy += potential;

            // Force scalar (without direction)
            double forceScalar =
                24.0 * m_epsilon *
                (2.0 * sigma12 * inv_r12 - sigma6 * inv_r6) * inv_r2;

            vec3 force = rij * forceScalar;

            atoms[i]->force += force;
            atoms[j]->force -= force;
        }
    }
}
