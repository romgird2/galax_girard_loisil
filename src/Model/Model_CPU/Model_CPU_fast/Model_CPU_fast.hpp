#ifdef GALAX_MODEL_CPU_FAST

#ifndef MODEL_CPU_FAST_HPP_
#define MODEL_CPU_FAST_HPP_

#include "../Model_CPU.hpp"
#include <cmath>
#include <iostream>



class Thread_Composition {
public:
    int clusters_ids[NB_CLUSTER_PER_THREAD];
    int complexity_total;

public:
    Thread_Composition() {
    }

    void init(int start)
    {
        complexity_total = 0;
        for(int i = 0;i != NB_CLUSTER_PER_THREAD;++i)
        {
            clusters_ids[i] = start+i;
        }
    }

    void calculate_properties(Particule *particules,Cluster *clusters)
    {
        for(int i = 0;i != NB_CLUSTER_PER_THREAD;++i)
        {
            int current_cluster_id = clusters_ids[i];
            Cluster &current_cluster = clusters[current_cluster_id];
            current_cluster.calculate_properties(particules);
        }
    }



    void execute(Particule *particules,Cluster *clusters,float distance_matrix_clusters[NB_TOTAL_CLUSTER][NB_TOTAL_CLUSTER])
    {
        complexity_total = 0;
        for(int i = 0;i != NB_CLUSTER_PER_THREAD;++i)
        {
            int current_cluster_id = clusters_ids[i];
            Cluster &current_cluster = clusters[current_cluster_id];
            current_cluster.complexity = current_cluster.nb_particules*current_cluster.nb_particules;

            // intra-cluster calculations
            for(int j = 0;j != current_cluster.nb_particules;++j)
            {
                int current_particule_id = current_cluster.particules[j];
                Particule &current_particule = particules[current_particule_id];
                current_particule.resetAcceleration();
                for(int k = 0;k != current_cluster.nb_particules;++k)
                {
                    if(k == j) continue;
                    int target_particule_id = current_cluster.particules[k];
                    Particule &target_particule = particules[target_particule_id];
                    Vector3 diff = target_particule.position-current_particule.position;
                    float distanceSquared = diff.normSquared();
                    if(distanceSquared < 1.0)
                        distanceSquared = 10.0;
                    else
                    {
                        distanceSquared = std::sqrt(distanceSquared);
                        distanceSquared = 10 /(distanceSquared*distanceSquared*distanceSquared);
                    }
                    current_particule.acceleration += diff*distanceSquared*target_particule.mass;

                }
            }

            // extra cluster calculations
            for(int target_cluster_id = 0;target_cluster_id != NB_TOTAL_CLUSTER;++target_cluster_id)
            {
                if(target_cluster_id == current_cluster_id) continue;
                Cluster &target_cluster = clusters[target_cluster_id];

                if(false)
                {

                }
                else
                {
                    current_cluster.complexity += current_cluster.nb_particules*target_cluster.nb_particules;
                    for(int j = 0;j != current_cluster.nb_particules;++j)
                    {
                        int current_particule_id = current_cluster.particules[j];
                        Particule &current_particule = particules[current_particule_id];
                        for(int k = 0;k != target_cluster.nb_particules;++k)
                        {
                            int target_particule_id = target_cluster.particules[k];
                            Particule &target_particule = particules[target_particule_id];
                            Vector3 diff = target_particule.position-current_particule.position;
                            float distanceSquared = diff.normSquared();
                            if(distanceSquared < 1.0)
                                distanceSquared = 2.0;
                            else
                            {
                                distanceSquared = std::sqrt(distanceSquared);
                                distanceSquared = 2 /(distanceSquared*distanceSquared*distanceSquared);
                            }
                            current_particule.acceleration += diff*distanceSquared*target_particule.mass;
                        }
                }
                }
            }
            complexity_total += current_cluster.complexity;
        }
    }
};



class Model_CPU_fast : public Model_CPU
{
public:
    Model_CPU_fast(const Initstate& initstate, Particles& particles);

    virtual ~Model_CPU_fast() = default;

    virtual void step();
    void balance_threads();

    Thread_Composition threads[NB_THREAD];
    Cluster clusters[NB_TOTAL_CLUSTER];
    Particule particules[NB_PARTICLES];
    int first_empty_cluster;

    float distance_matrix_clusters[NB_TOTAL_CLUSTER][NB_TOTAL_CLUSTER];

};






#endif // MODEL_CPU_FAST_HPP_

#endif // GALAX_MODEL_CPU_FAST
