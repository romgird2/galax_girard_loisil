

#ifndef MODEL_CPU_FAST_HPP_
#define MODEL_CPU_FAST_HPP_

#include "../Model_CPU.hpp"
#include <cmath>
#include <iostream>





class Model_CPU_fast : public Model_CPU
{
public:
    Model_CPU_fast(const Initstate& initstate, Particule *particles);

    virtual ~Model_CPU_fast() = default;

    virtual void step();

    Vector3 origin;
    int ***particules_repartition;
    int **nb_particules_repartition;

    OctTree root;

};






#endif // MODEL_CPU_FAST_HPP_

