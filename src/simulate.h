/** \file simulate.h
 *
 * Functions for simulating phylogenies on contact networks.
 */

#ifndef SIMULATE_H
#define SIMULATE_H

#include <igraph/igraph.h>
#include <gsl/gsl_rng.h>

/** Simulate a phylogenetic tree from a contact network.
 *
 * This uses a simple algorithm. Let R be the sum of the rates of all
 * sero-discordant edges in the network, that is, edges where one endpoint is
 * infected and the other is susceptible. Whenever a new node is infected, we
 * update R accordingly.
 *
 * The epidemic is started by infecting a random node. The time to the next
 * infection is chosen from Exponential(R), a sero-discordant edge is chosen at
 * random, and the susceptible node along it is infected. The phylogeny is
 * constructed by keeping track of who transmitted to who.
 *
 * If the contact network is not connected, the phylogeny will only be
 * simulated on a (random) connected component. The stop_time and stop_tips
 * parameters are maximums, so the returned tree may have fewer nodes and/or be
 * shorter. This function produces a complete phylogeny (no tip sampling;
 * removed tips are included at their removal time).
 *
 * The contact network should have the following edge attributes:
 * * "transmit": the transmission rate along the edge
 *
 * The contact network should have the following vertex attributes:
 * * "remove": the removal rate of the node once it becomes infected (if the
 * node can recover or die, the removal rate is the sum of the recovery and
 * death rates)
 * 
 * If a non-positive value is passed for either the stop_time or stop_nodes
 * parameters, no limit will be imposed. So if stop_nodes <= 0, the simulation
 * will continue until no more sero-discordant edges exist in the network. This
 * does not necessarily mean every node will be infected, if the network is not
 * connected.
 *
 * \param[in] net the contact network, with rates defined on each edge
 * \param[in] rng the GSL random generator object
 * \param[in] stop_time maximum amount of time to run the simulation for, <= 0 means no limit
 * \param[in] stop_nodes maximum number of nodes to infect, <= 0 means no limit
 * \param[in] numeric_ids 1 if the "id" attribute of the network is numeric
 * \param[in] tree an uninitialized igraph_t object
 * \return a phylogeny simulated over the contact network
 */
void simulate_phylogeny(igraph_t *tree, igraph_t *net, gsl_rng *rng, double
        stop_time, int stop_nodes, int numeric_ids);

#endif
