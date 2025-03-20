#include "Model.hpp"
#include <cassert>
#include <numeric>
#include <cmath>
#include <iostream>
#include <algorithm>

Model
::Model(const Initstate& initstate, Particule *particules)
: initstate(initstate),
  particules(particules)
{
}

std::tuple<float, float, float> Model
::compareParticlesState(const Model& referenceModel, bool returnRelativeDistances)
{
    // Compute the average distance between the particles in the two datasets.
    // We could also do a relative error, but given that we expect extremely
    // close results regardless of the Model used, this should be good enough.
    std::vector<float> distances(NB_PARTICLES, 0.0); // Maybe make this static
    std::vector<float> relative_distances(NB_PARTICLES, 0.0); // Maybe make this static

    // Should we parallelize?
    for(size_t p = 0; p < NB_PARTICLES; ++p)
    {
        distances[p] = std::sqrt((particules[p].position.x - referenceModel.particules[p].position.x) * (particules[p].position.x - referenceModel.particules[p].position.x) +
                                 (particules[p].position.y - referenceModel.particules[p].position.y) * (particules[p].position.y - referenceModel.particules[p].position.y) +
                                 (particules[p].position.z - referenceModel.particules[p].position.z) * (particules[p].position.z - referenceModel.particules[p].position.z));

        relative_distances[p] = distances[p] * std::sqrt((referenceModel.particules[p].position.x) * (referenceModel.particules[p].position.x) +
                                                         (referenceModel.particules[p].position.y) * (referenceModel.particules[p].position.y) +
                                                         (referenceModel.particules[p].position.z) * (referenceModel.particules[p].position.z));
    }

    if(returnRelativeDistances)
    {
        auto minmax = minmax_element(distances.cbegin(), distances.cend());
        return {*minmax.first, *minmax.second, std::accumulate(distances.cbegin(), distances.cend(), 0.0) / distances.size()};
    }
    else
    {
        auto minmax = minmax_element(relative_distances.cbegin(), relative_distances.cend());
        return {*minmax.first, *minmax.second, std::accumulate(relative_distances.cbegin(), relative_distances.cend(), 0.0) / relative_distances.size()};
    }
}

