#include <cmath>

#include "Model_CPU.hpp"

Model_CPU::Model_CPU(const Initstate& initstate, Particule *particules)
: Model(initstate, particules)
{
    for (int i = 0; i < NB_PARTICLES; i++)
	{
        particules[i].position.set(initstate.positionsx[i],initstate.positionsy[i],initstate.positionsz[i]);
        particules[i].velocity.set(initstate.velocitiesx[i]*0.2,initstate.velocitiesy[i]*0.2,initstate.velocitiesz[i]*0.2);
        particules[i].mass = initstate.masses[i];
    }
}
