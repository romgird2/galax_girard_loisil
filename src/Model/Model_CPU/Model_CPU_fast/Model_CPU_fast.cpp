#ifdef GALAX_MODEL_CPU_FAST

#include <cmath>

#include "Model_CPU_fast.hpp"

#include <xsimd/xsimd.hpp>
#include <omp.h>
#include <ranges>

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

bool particules_x_sorting (int i,int j) { return (i<j); }

bool particules_y_sorting (int i,int j) { return (i<j); }

bool particules_z_sorting (int i,int j) { return (i<j); }

void separate_particles(const Particule* particles,Cluster *clusters,std::vector<int> &particles_id,std::vector<int> &clusters_id,int cord_index)
{
    if(clusters_id.size() == 0) return;
    if(clusters_id.size() == 1)
    {
        if(particles_id.size() > MAX_PARTICULES_PER_CLUSTER)
            throw std::exception();
        Cluster &cluster = clusters[clusters_id.at(0)];
        cluster.nb_particules = particles_id.size();
        for(int i = 0;i != particles_id.size();++i)
            cluster.particules[i] = particles_id[i];

        return;
    }

    auto f = [particles] (int i,int j, auto field) {
            return particles[i].position.*field < particles[j].position.*field;
        };

    switch(cord_index)
    {
    case 0:
        std::sort(particles_id.begin(),particles_id.end(),[particles,f ](int i, int j) {return f(i,j,&Vector3::x);});
        break;
    case 1:
        std::sort(particles_id.begin(),particles_id.end(),[particles,f ](int i, int j) {return f(i,j,&Vector3::y);});
        break;
    default:
        std::sort(particles_id.begin(),particles_id.end(),[particles,f ](int i, int j) {return f(i,j,&Vector3::z);});
        break;
    }



    std::size_t const half_size = particles_id.size() / 2;
    std::vector<int> split_lo(particles_id.begin(), particles_id.begin() + half_size);
    std::vector<int> split_hi(particles_id.begin() + half_size, particles_id.end());

    std::size_t const half_size_clusters = clusters_id.size() / 2;
    std::vector<int> split_lo_clusters(clusters_id.begin(), clusters_id.begin() + half_size_clusters);
    std::vector<int> split_hi_clusters(clusters_id.begin() + half_size_clusters, clusters_id.end());

    separate_particles(particles,clusters,split_lo,split_hi_clusters,cord_index+1%3);
    separate_particles(particles,clusters,split_hi,split_lo_clusters,cord_index+1%3);
}


Model_CPU_fast::Model_CPU_fast(const Initstate& initstate, Particles& particles)
: Model_CPU(initstate, particles)
{
    for(int i = 0;i != NB_PARTICLES;++i)
    {
        Particule &particule = particules[i];
        particule.position.set(initstate.positionsx.at(i),initstate.positionsy.at(i),initstate.positionsz.at(i));
        particule.velocity.set(initstate.velocitiesx.at(i)*0.1f,initstate.velocitiesy.at(i)*0.1f,initstate.velocitiesz.at(i)*0.1f);
        particule.mass = initstate.masses.at(i);
    }



    for(int i = 0;i != NB_TOTAL_CLUSTER;++i)
        threads[i].init(NB_CLUSTER_PER_THREAD*i);


    std::vector<int> particules_id(NB_PARTICLES),clusters_id(NB_TOTAL_CLUSTER);

    for(int i = 0;i != NB_PARTICLES;++i)
        particules_id[i] = i;

    for(int i = 0;i != NB_TOTAL_CLUSTER;++i)
        clusters_id[i] = i;

    separate_particles(particules,clusters,particules_id,clusters_id,0);

    for(int i = 0;i != NB_TOTAL_CLUSTER;++i)
    {
        Cluster &cluster = clusters[i];
        cluster.init_total_mass(particules);
        distance_matrix_clusters[i][i] = 0;
    }
}

void Model_CPU_fast::balance_threads()
{
    int min_complexity = threads->complexity_total,min_complexity_id = 0,
        max_complexity = threads->complexity_total,max_complexity_id = 0;

    for(int i = 1;i != NB_THREAD;++i)
    {
        int &complexity = threads[i].complexity_total;
        if(complexity < min_complexity)
        {
            min_complexity = complexity;
            min_complexity_id = i;
        }
        if(complexity > max_complexity)
        {
            max_complexity = complexity;
            max_complexity_id = i;
        }
    }

    Thread_Composition &thread_not_enough = threads[min_complexity_id];
    Thread_Composition &thread_too_much = threads[max_complexity_id];
    int cluster_id_id_min = 0,
        cluster_id_id_max = 0;
    int cluster_complexity_min = clusters[thread_too_much.clusters_ids[cluster_id_id_min]].complexity;
    int cluster_complexity_max = clusters[thread_not_enough.clusters_ids[cluster_id_id_max]].complexity;

    for(int i = 1;i != NB_CLUSTER_PER_THREAD;++i)
    {

        int &cluster_id_not_enough = thread_not_enough.clusters_ids[i];
        Cluster &cluster_not_enough = clusters[cluster_id_not_enough];

        int &cluster_id_too_much = thread_too_much.clusters_ids[i];
        Cluster &cluster_too_much = clusters[cluster_id_too_much];

        if(cluster_too_much.complexity > cluster_complexity_max)
        {
            cluster_complexity_max = cluster_too_much.complexity;
            cluster_id_id_max = i;
        }
        if(cluster_not_enough.complexity < cluster_complexity_min)
        {
            cluster_complexity_min = cluster_not_enough.complexity;
            cluster_id_id_min = i;
        }
    }
    std::swap(thread_too_much.clusters_ids[cluster_id_id_max],thread_not_enough.clusters_ids[cluster_id_id_min]);

}

bool balance_between(Cluster &current_cluster,Cluster &target_cluster,Particule *particules)
{
    if(current_cluster.nb_particules > MIN_PARTICULES_PER_CLUSTER && target_cluster.nb_particules+1 < MAX_PARTICULES_PER_CLUSTER)
    {
        for(int i = 0;i != current_cluster.nb_particules;++i)
        {
            Particule &current_particule = particules[current_cluster.particules[i]];
            float delta_x = std::abs(target_cluster.center.x-current_particule.position.x);
            float delta_y = std::abs(target_cluster.center.y-current_particule.position.y);
            float delta_z = std::abs(target_cluster.center.z-current_particule.position.z);
            float distance_from_target = std::max(delta_x,std::max(delta_y,delta_z));
            if(distance_from_target+MIN_DISTANCE_FUSE < current_cluster.distance_from_center[i])
            {
                target_cluster.steal(current_cluster,i,particules);
                return true;
            }
        }
    }
    return false;
}

void Model_CPU_fast
::step()
{
    #pragma omp parallel for
    for(int i = 0;i != NB_THREAD;++i)
    {
        Thread_Composition& thread = threads[i];
        thread.calculate_properties(particules,clusters);
    }

    #pragma omp parallel for
    for(int i = 0;i != NB_TOTAL_CLUSTER;++i)
    {
        Cluster &current_cluster = clusters[i];
        for(int j = 0;j != i;++j)
        {
            Cluster &target_cluster = clusters[j];
            float distance = (current_cluster.center-target_cluster.center).norm();
            distance_matrix_clusters[i][j] = distance;
            distance_matrix_clusters[j][i] = distance;
        }
    }

    bool updated[NB_TOTAL_CLUSTER];
    for(int i = 0;i != NB_TOTAL_CLUSTER;++i)
        updated[i] = true;
    for(int i = 0;i != NB_TOTAL_CLUSTER;++i)
    {
        Cluster &current_cluster = clusters[i];
        for(int target_cluster_id = 0;target_cluster_id != i;++target_cluster_id)
        {
            if(updated[i] && updated[target_cluster_id])
            {
                Cluster &target_cluster = clusters[target_cluster_id];
                if(distance_matrix_clusters[i][target_cluster_id] < std::max(current_cluster.radius,target_cluster.radius))
                {
                    if(current_cluster.nb_particules > target_cluster.nb_particules)
                    {
                        if(balance_between(current_cluster,target_cluster,particules))
                        {
                            float distance = (target_cluster.center-current_cluster.center).norm();
                            distance_matrix_clusters[i][target_cluster_id] = distance;
                            distance_matrix_clusters[target_cluster_id][i] = distance;
                            updated[i] = false;
                            goto stop_updating_cluster;
                        }
                        if(balance_between(target_cluster,current_cluster,particules))
                        {
                            float distance = (target_cluster.center-current_cluster.center).norm();
                            distance_matrix_clusters[i][target_cluster_id] = distance;
                            distance_matrix_clusters[target_cluster_id][i] = distance;
                            updated[target_cluster_id] = false;
                            goto stop_updating_cluster;
                        }
                    }
                    else
                    {
                        if(balance_between(target_cluster,current_cluster,particules))
                        {
                            float distance = (target_cluster.center-current_cluster.center).norm();
                            distance_matrix_clusters[i][target_cluster_id] = distance;
                            distance_matrix_clusters[target_cluster_id][i] = distance;
                            updated[target_cluster_id] = false;
                            goto stop_updating_cluster;
                        }
                        if(balance_between(current_cluster,target_cluster,particules))
                        {
                            float distance = (target_cluster.center-current_cluster.center).norm();
                            distance_matrix_clusters[i][target_cluster_id] = distance;
                            distance_matrix_clusters[target_cluster_id][i] = distance;
                            updated[i] = false;
                            goto stop_updating_cluster;
                        }
                    }
                }
            }
        }
        stop_updating_cluster:;
    }


    #pragma omp parallel for
    for(int i = 0;i < NB_THREAD;++i)
    {
        Thread_Composition& thread = threads[i];
        thread.execute(particules,clusters,distance_matrix_clusters);
    }



    balance_threads();

    float max_radius = 0;
    for(int i = 0;i != NB_TOTAL_CLUSTER;++i)
    {
        if(clusters[i].radius > max_radius)
            max_radius = clusters[i].radius;
    }
std::cout << "max radius is " << max_radius << std::endl;

    for(int i = 0;i != NB_PARTICLES;++i)
    {
        Particule& particule = particules[i];
        particule.velocity += particule.acceleration;
        particule.position += particule.velocity;
        particles.x.at(i) = particule.position.x;
        particles.y.at(i) = particule.position.y;
        particles.z.at(i) = particule.position.z;
    }
}

#endif // GALAX_MODEL_CPU_FAST
