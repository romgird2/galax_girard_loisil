#ifdef GALAX_MODEL_CPU_FAST

#include <cmath>

#include "Model_CPU_fast.hpp"

#include <xsimd/xsimd.hpp>
#include <omp.h>
#include <ranges>


Model_CPU_fast::Model_CPU_fast(const Initstate& initstate, Particule *particles)
: Model_CPU(initstate, particles)
{
    omp_set_num_threads(NB_THREADS);

    int sorted_x_indexes[NB_PARTICLES],sorted_y_indexes[NB_PARTICLES],sorted_z_indexes[NB_PARTICLES];
    for(int i = 0;i != NB_PARTICLES;++i)
    {
        sorted_x_indexes[i] = NB_PARTICLES-i-1;
        sorted_y_indexes[i] = i;
        sorted_z_indexes[i] = i;
    }
    std::sort(&sorted_x_indexes[0],&sorted_x_indexes[NB_PARTICLES],[this] (int i,int j) {return particules[i].position.x < particules[j].position.x;});
    std::sort(&sorted_y_indexes[0],&sorted_y_indexes[NB_PARTICLES],[this] (int i,int j) {return particules[i].position.y < particules[j].position.y;});
    std::sort(&sorted_z_indexes[0],&sorted_z_indexes[NB_PARTICLES],[this] (int i,int j) {return particules[i].position.z < particules[j].position.z;});
    origin = Vector3(particules[sorted_x_indexes[NB_PARTICLES/2]].position.x,
                    particules[sorted_y_indexes[NB_PARTICLES/2]].position.y,
                    particules[sorted_z_indexes[NB_PARTICLES/2]].position.z);

    OctTree::particules = particles;


    particules_repartition = new int**[NB_THREADS];
    for(int i = 0;i != NB_THREADS;++i)
    {
        particules_repartition[i] = new int*[NB_THREADS];
        for(int j = 0;j != NB_THREADS;++j)
            particules_repartition[i][j] = new int[NB_PARTICLES/NB_THREADS+NB_THREADS];
    }

    nb_particules_repartition = new int*[NB_THREADS];
    for(int i = 0;i != NB_THREADS;++i)
        nb_particules_repartition[i] = new int[NB_THREADS];

    root.init(Vector3(0,0,0),1);
}




void Model_CPU_fast
::step()
{

    root.clear();

    int *max_x = new int[NB_THREADS],*max_y = new int[NB_THREADS],*max_z = new int[NB_THREADS];
    int *min_x = new int[NB_THREADS],*min_y = new int[NB_THREADS],*min_z = new int[NB_THREADS];


    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();

        std::cout << "repartition initial " << thread_id << std::endl;

        int start = thread_id * (NB_PARTICLES / NB_THREADS);
        int end = (thread_id == NB_THREADS - 1) ? NB_PARTICLES : (thread_id + 1) * (NB_PARTICLES / NB_THREADS);

        int **particules_repartition_thread = particules_repartition[thread_id];
        int *nb_particules_repartition_thread = nb_particules_repartition[thread_id];
        for(int i = 0;i != 8;++i)
            nb_particules_repartition_thread[i] = 0;
        int &max_x_thread = max_x[thread_id],&max_y_thread = max_y[thread_id],&max_z_thread = max_z[thread_id];
        int &min_x_thread = min_x[thread_id],&min_y_thread = min_y[thread_id],&min_z_thread = min_z[thread_id];

        max_x_thread = 0;max_y_thread = 0;max_z_thread = 0;
        min_x_thread = 0;min_y_thread = 0;min_z_thread = 0;
        for(int i = start;i != end;++i)
        {
            Vector3& position = particules[i].position;
            if(position.x > max_x_thread) max_x_thread = position.x;
            if(position.x < min_x_thread) min_x_thread = position.x;

            if(position.y > max_y_thread) max_y_thread = position.y;
            if(position.y < min_y_thread) min_y_thread = position.y;

            if(position.z > max_z_thread) max_z_thread = position.z;
            if(position.z < min_z_thread) min_z_thread = position.z;

            int part_index = (position.x > origin.x) | ((position.y > origin.y)<<1) | ((position.z > origin.z)<<2);
            particules_repartition_thread[part_index][nb_particules_repartition_thread[part_index]++] = i;
        }
    }

    int global_max_x = 0, global_max_y = 0, global_max_z = 0;
    int global_min_x = 0, global_min_y = 0, global_min_z = 0;

    for (int i = 0; i < NB_THREADS; ++i) {
        if (max_x[i] > global_max_x) global_max_x = max_x[i];
        if (max_y[i] > global_max_y) global_max_y = max_y[i];
        if (max_z[i] > global_max_z) global_max_z = max_z[i];

        if (min_x[i] < global_min_x) global_min_x = min_x[i];
        if (min_y[i] < global_min_y) global_min_y = min_y[i];
        if (min_z[i] < global_min_z) global_min_z = min_z[i];
    }

    int max_abs_diff = std::max({
        global_max_x - origin.x,
        global_max_y - origin.y,
        global_max_z - origin.z,
        origin.x - global_min_x,
        origin.y - global_min_y,
        origin.z - global_min_z
    });

    root.init(origin,max_abs_diff);
    root.unleaf();


    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        std::cout << "repartition into 8 " << thread_id << std::endl;
        OctTree &target = root.children[thread_id];

        for(int i = 0;i != NB_THREADS;++i)
        {
            int *particules_repartition_thread = particules_repartition[i][thread_id];
            int nb_particule_repartition_thread = nb_particules_repartition[i][thread_id];
            for(int j = 0;j != nb_particule_repartition_thread;++j)
            {
                int particule_index = particules_repartition_thread[j];
                target.insert(particules[particule_index].position,particule_index);
            }
        }
        target.pre_compute();
    }

    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        std::cout << "calculate 8 " << thread_id << std::endl;

        OctTree &target = root.children[thread_id];
        target.compute_acceleration();
    }

    for(int i = 0;i != NB_PARTICLES;++i)
    {
        Particule& particule = particules[i];
        particule.velocity += particule.acceleration;
        //particule.position += particule.velocity;
    }
}

#endif // GALAX_MODEL_CPU_FAST
