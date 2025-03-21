#include "Particles.hpp"

Particule *OctTree::particules;



Particles::
Particles(const int n_particles)
: x(n_particles),
  y(n_particles),
  z(n_particles)
{
}
