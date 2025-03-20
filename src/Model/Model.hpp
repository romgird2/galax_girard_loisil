#ifndef MODEL_HPP_
#define MODEL_HPP_

#include <vector>
#include <tuple>

#include "../Initstate.hpp"
#include "../Particles.hpp"

class Model
{
public:
    const Initstate& initstate;

    Particule *particules;

public:
    Model(const Initstate& initstate, Particule *particules);

    std::tuple<float, float, float> compareParticlesState(const Model& reference, bool returnRelativeDistances = false);

    virtual ~Model() = default;

    virtual void step() = 0;
};

#endif // MODEL_HPP
