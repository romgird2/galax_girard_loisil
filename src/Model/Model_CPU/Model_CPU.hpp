#ifndef MODEL_CPU_HPP_
#define MODEL_CPU_HPP_

#include "../Model.hpp"

class Model_CPU : public Model
{

public:
    Model_CPU(const Initstate& initstate, Particule *particules);

    virtual ~Model_CPU() = default;

    virtual void step() = 0;
};
#endif // MODEL_CPU_NAIVE_HPP_
