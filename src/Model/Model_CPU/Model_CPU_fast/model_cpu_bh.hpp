#ifndef MODEL_CPU_BH_H
#define MODEL_CPU_BH_H

#include "../Model_CPU.hpp"
#include "../../Body.hpp"
#include <vector>
#include <array>
#include "../BHTree.hpp"

/**
 * @brief The Model_CPU_BH class
 * Barnes-Hut NBody simulation
 */
class Model_CPU_BH : public Model_CPU
{
public:
    Model_CPU_BH(const Initstate& initstate, Particles& particles);

    virtual ~Model_CPU_BH() = default;

    virtual void step();
private:
    std::vector<Body> bodies;
    // Repartition of particles for each thread for tree creation
    std::array<std::vector<Body>,8> bodies_per_thread;
    // Tree created by each thread
    std::array<std::unique_ptr<BHTree>,8> trees;
    double radius;
    Vector3 last_mass_center = {0.0, 0.0, 0.0};
};

#endif // MODEL_CPU_BH_H
