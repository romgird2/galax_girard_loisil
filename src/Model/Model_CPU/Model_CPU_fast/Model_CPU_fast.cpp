#ifdef GALAX_MODEL_CPU_FAST

#include <cmath>

#include "Model_CPU_fast.hpp"

#include <xsimd/xsimd.hpp>
#include <omp.h>

namespace xs = xsimd;
using b_type = xs::batch<float, xs::avx2>;
void computeBounds(const Particule* particles, int num_particles, Vector3& min_pos, Vector3& max_pos) {
    if (num_particles == 0) return;

    min_pos = particles[0].position;
    max_pos = particles[0].position;

    for (int i = 1; i < num_particles; ++i) {
        const Vector3& pos = particles[i].position;
        min_pos.x = std::min(min_pos.x, pos.x);
        min_pos.y = std::min(min_pos.y, pos.y);
        min_pos.z = std::min(min_pos.z, pos.z);
        max_pos.x = std::max(max_pos.x, pos.x);
        max_pos.y = std::max(max_pos.y, pos.y);
        max_pos.z = std::max(max_pos.z, pos.z);
    }
}


Model_CPU_fast::Model_CPU_fast(const Initstate& initstate, Particles& particles)
: Model_CPU(initstate, particles)
{
    for(int i = 0;i != NB_PARTICLES;++i)
    {
        Particule &particule = particules[i];
        particule.position.set(initstate.positionsx.at(i),initstate.positionsy.at(i),initstate.positionsz.at(i));
        particule.velocity.set(initstate.velocitiesx.at(i),initstate.velocitiesy.at(i),initstate.velocitiesz.at(i));
        particule.mass = initstate.masses.at(i);
    }

    for(int i = 0;i != NB_MAX_CLUSTER;++i)
    {
        Cluster &cluster = clusters[i];

    }

    for(int i = 0;i != NB_MAX_CLUSTER;++i)
        threads[i].init(NB_MAX_CLUSTER_PER_THREAD*i);

    Vector3 minBound,maxBound;
    computeBounds(particules,NB_PARTICLES,minBound,maxBound);

    int nb_regions = std::floor(std::cbrt(NB_MAX_CLUSTER));

    float invdx = nb_regions/(maxBound.x - minBound.x);
    float invdy = nb_regions/(maxBound.y - minBound.y);
    float invdz = nb_regions/(maxBound.z - minBound.z);

    for (int i = 0; i != NB_PARTICLES; ++i) {
        const Vector3& pos = particules[i].position;

        int ix = static_cast<int>((pos.x - minBound.x) *invdx);
        ix = std::clamp(ix, 0, nb_regions - 1);

        int iy = static_cast<int>((pos.y - minBound.y) *invdy);
        iy = std::clamp(iy, 0, nb_regions - 1);

        int iz = static_cast<int>((pos.z - minBound.z) *invdz);
        iz = std::clamp(iz, 0, nb_regions - 1);

        int cluster_id = iz * (nb_regions * nb_regions) + iy * nb_regions + ix;


        Cluster* cluster = &clusters[cluster_id];
        while(cluster->nb_particules == MAX_PARTICULES_PER_CLUSTER)
        {
            cluster_id = (cluster_id+1)%NB_MAX_CLUSTER;
            cluster = &clusters[cluster_id];
        }
        cluster->particules[cluster->nb_particules++] = i;
    }
}

void Model_CPU_fast
::step()
{
    #pragma omp parallel for
    for(int i = 0;i < NB_THREAD;++i)
    {
        Thread_Composition* thread = &threads[i];
        thread->execute(particules,clusters);
    }


    for(int i = 0;i != NB_PARTICLES;++i)
    {
        Particule& particule = particules[i];
        particule.velocity += particule.acceleration * 2;
        particule.position += particule.velocity * 0.1;
        particles.x.at(i) = particule.position.x;
        particles.y.at(i) = particule.position.y;
        particles.z.at(i) = particule.position.z;
    }
}

#endif // GALAX_MODEL_CPU_FAST
