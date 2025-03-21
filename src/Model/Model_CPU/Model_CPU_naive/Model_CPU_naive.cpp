#include <cmath>

#include "Model_CPU_naive.hpp"

Model_CPU_naive
::Model_CPU_naive(const Initstate& initstate, Particule *particles)
: Model_CPU(initstate, particles)
{
}

void Model_CPU_naive
::step()
{
    for(int i = 0;i != NB_PARTICLES;++i)
    {
        particules[i].acceleration.set(0,0,0);
    }

    for (int i = 0; i < NB_PARTICLES; i++)
	{
        for (int j = 0; j < NB_PARTICLES; j++)
		{
			if(i != j)
			{
                Vector3 diff = particules[i].position-particules[j].position;
                float dij = diff.normSquared();

				if (dij < 1.0)
				{
                    dij = 2.0;
				}
				else
				{
					dij = std::sqrt(dij);
                    dij = 2.0 / (dij * dij * dij);
				}

                particules[i].acceleration += diff * dij * particules[j].mass;
			}
		}
	}

    for (int i = 0; i < NB_PARTICLES; i++)
	{
        particules[i].velocity += particules[i].acceleration;
        particules[i].position += particules[i].velocity;
	}
}
