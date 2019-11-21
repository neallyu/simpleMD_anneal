#ifndef NEIGHBORLIST_H
#define NEIGHBORLIST_H

#include <vector>
#include "particle.hpp"
#include "ensemble_reduced_unit.hpp"

using namespace std;

class Neighborlist {

friend class Ensemble;

public:
    Neighborlist(vector<Particle> &, double box, double _rlist2);

    void update_neighbor_list(vector<Particle> &);

private:
    double calc_distance2(Particle& particle1, Particle& particle2);

    vector<vector <int> > nlist;
    const double BOX;
    const double rlist2;
};


Neighborlist::Neighborlist(vector<Particle> &ensemble, double box, double _rlist2): 
    BOX(box), rlist2(_rlist2) {
        nlist.resize(ensemble.size());
    }

void Neighborlist::update_neighbor_list(vector<Particle> &ensemble) {
    for (int i = 0; i < nlist.size(); i++)
    {
        nlist[i].resize(0);
    }

// Atoms are not double counted in the neighbor list. That is, when atom j
// is on atom i's list, the opposite is not true.
    for (int i = 0; i < nlist.size() - 1; ++i) {
        for (int j = i + 1; j < nlist.size(); ++j) {
            if (calc_distance2(ensemble[i], ensemble[j]) < rlist2) {
                nlist[i].push_back(j);
            }
        }
    }
}


double Neighborlist::calc_distance2(Particle& particle1, Particle& particle2) {
    double dx = particle1.pos_x - particle2.pos_x;
    double dy = particle1.pos_y - particle2.pos_y;
    double dz = particle1.pos_z - particle2.pos_z;

    // periodic boundary conditions
    dx -= BOX * round(dx / BOX);
    dy -= BOX * round(dy / BOX);
    dz -= BOX * round(dz / BOX);

    return dx * dx + dy * dy + dz * dz;
}


#endif