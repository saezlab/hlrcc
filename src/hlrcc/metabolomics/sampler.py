#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import numpy as np
import random as rand
from pandas import DataFrame
from scipy.linalg.lapack import dgesdd
from scipy.linalg.blas import dgemm
from framed import FBA, FVA


def set_null_space(matrix, tol=1e-12):
    # Single value decomposition
    u, s, vt, info = dgesdd(matrix)

    # Singular values samller than tol are considered zero
    rank = (s > tol).sum()

    # Null space are the columns of v corresponding to the zero singular values
    n = vt[rank:].transpose()

    return n, n.transpose()


def get_warmup_points(model, reactions):
    # Calculate warm-up points running FVA
    fva, simulations = FVA(model, return_simulations=True)

    # Get upper and lower bounds
    lb, ub = np.array(zip(*[fva[r] for r in reactions]))

    # Warm-up points
    warmup_points = simulations.loc[reactions].values

    return lb, ub, warmup_points


def project_onto_null_space(a, n, nt, lb, ub, n_steps_projection, reactions, tol=1e-4, verbose=0):
    # Do matrix multuplication
    a = dgemm(1.0, n, dgemm(1.0, nt, a)).transpose()[0]

    i, accept = 0, True
    while any(a < lb - tol) or any(a > ub + tol):
        print '[WARNING] Bounds reseted: %d times' % (i+1)
        print '\t lower bound: ', reactions[a < lb - tol]
        print '\t upper bound: ', reactions[a > ub + tol]

        # Fix bounds violations
        lb_violations = a < lb - tol
        if any(lb_violations):
            a[lb_violations] = lb[lb_violations]

        ub_violations = a > ub
        if any(ub_violations):
            a[ub_violations] = ub[ub_violations]

        # Do matrix multiplication
        a = dgemm(1.0, n, dgemm(1.0, nt, a)).transpose()[0]

        # Verbose
        if i > 10:
            accept, n_steps_projection = False, max(25, n_steps_projection - 100)
            # Verbose
            if verbose > 1:
                print '[WARNING] Point discarded! lb: %r, ub: %r' % (any(lb_violations), any(ub_violations))
            break

        i += 1

    if i == 0:
        n_steps_projection = min(25, n_steps_projection + 25)

    return accept, n_steps_projection, a.copy()


def get_point(current_point, centre_point, samples, lb, ub, tol=1e-7, verbose=0):
    accept = True

    # Get random point index from the samples
    r_index = rand.randint(0, len(samples)-1)

    # Calculate direction from the centre point to the sampled point
    u = samples[r_index] - centre_point
    u /= np.linalg.norm(u, 2)

    # Check move forward or backward
    neg_dir, pos_dir = u < -tol, u > tol
    movable = np.bitwise_or(neg_dir, pos_dir)

    # If movable is all False continue cycle
    if np.all(np.logical_not(movable)):
        print '[WARNING] Not movable: ', movable
        accept = False

    else:
        # Get distance to the lower and upper bound for the movable direactions
        min_step_vec = np.concatenate([
            -(current_point[pos_dir] - lb[pos_dir]) / u[pos_dir],
            (ub[neg_dir] - current_point[neg_dir]) / u[neg_dir]
        ])

        max_step_vec = np.concatenate([
            (ub[pos_dir] - current_point[pos_dir]) / u[pos_dir],
            -(current_point[neg_dir] - lb[neg_dir]) / u[neg_dir]
        ])

        min_step = min(np.max(min_step_vec), 0.0)
        max_step = max(np.min(max_step_vec), 0.0)

        r_step = rand.random()
        step_size = min_step + r_step * (max_step - min_step)

        # Make the step
        current_point[movable] += step_size * u[movable]

        # Verbose
        if verbose > 1:
            print '[INFO] #Samples: %d, Step size = %.5f, Random index: %d, Random step: %.5f' % (len(samples) - len(current_point) * 2, step_size, r_index, r_step)

    return accept, current_point


def sample(model, n_samples=1000, n_steps=50, n_steps_projection=25, verbose=0):
    # Set general variables
    in_null_space, n_pts_discarded = False, 0

    # Get reactions names
    reactions = np.array(model.reactions.keys())

    # Get stoichiometic matrix as numpy arrays
    S = model.get_stoichiometric_matrix().loc[:, reactions].values

    # Compute the null space of the model
    n, nt = set_null_space(S)

    # Calculate warm-up points
    lb, ub, warmup_points = get_warmup_points(model, reactions)

    # Calculate centre point
    centre_point = warmup_points.mean(1)

    # Move warm-up point towards the centre
    warmup_points = np.multiply(warmup_points, 0.33) + np.multiply(np.resize(centre_point, (warmup_points.shape[0], 1)), 0.67).dot(np.ones((1, warmup_points.shape[1])))

    # Initialise samples with warm-up points
    samples = warmup_points.transpose().copy()
    while len(samples) < n_samples + warmup_points.shape[1]:
        # Check discarded samples
        if n_pts_discarded >= 1000:
            print '[ERROR] Too many points discarded, sampling will be aborted: %d' % n_pts_discarded
            return None

        elif n_pts_discarded >= 100:
            print '[WARNING] Discarded point: %d' % n_pts_discarded

        # Pick the current point
        if not in_null_space:
            r_cp_index = rand.randint(0, len(samples)-1)
            current_point = samples[r_cp_index].copy()
            in_null_space = True

        step = 1
        while step <= n_steps:
            # Get next point
            accept, current_point = get_point(current_point, centre_point, samples, lb, ub, verbose=verbose)
            if not accept:
                print '[ERROR] get_point: New point discarded: ', current_point
                n_pts_discarded += 1
                step, in_null_space = n_steps, False
                continue

            # Correct floating point errors by projecting onto N
            if step % n_steps_projection == 0 and in_null_space:
                accept, n_steps_projection, current_point = project_onto_null_space(current_point, n, nt, lb, ub, n_steps_projection, reactions, verbose=verbose)
                if not accept:
                    print '[WARNING] project_onto_null_space: New point discarded: ', current_point
                    n_pts_discarded += 1
                    in_null_space = False
                    continue

            step += 1

        if not in_null_space:
            continue

        # Final projection
        if step % n_steps_projection != 0:
            accept, n_steps_projection, current_point = project_onto_null_space(current_point, n, nt, lb, ub, n_steps_projection, reactions, verbose=verbose)

            if not accept or not in_null_space:
                print '[WARNING] Final projection: Discarded point: %d' % n_pts_discarded
                n_pts_discarded += 1
                in_null_space = False
                continue

        # Update centre point
        if len(samples) % 25 == 0:
            accept, n_steps_projection, centre_point = project_onto_null_space(centre_point, n, nt, lb, ub, n_steps_projection, reactions, verbose=verbose)

        centre_point = (centre_point * len(samples) + current_point) / (len(samples) + 1)

        # Store current point
        samples = np.vstack((samples, current_point.copy()))

        # Verbose
        if verbose > 0:
            print '[INFO] Samples: %d' % (len(samples) - warmup_points.shape[1])

    if verbose > 0:
        print '[INFO] Sampling done!: %d points discarded' % n_pts_discarded

    return DataFrame(samples[warmup_points.shape[1]:], columns=reactions)


def fix_futile_cycles(model, samples, verbose=1):
    """ Fixes the distribution values of the reactions in futile cycles, by fixing all the other reactions values
        and minimising the futile cycle reactions values.

    Arguements:
        model: Model -- a metabolic constrain-based model
        samples: DataFrame -- pandas DataFrame returned by the sample function - rows samples, columns reactions

    Returns:
        solution: DataFrame -- pandas DataFrame same structure as the samples argument
    """
    m = model.deepcopy()

    # Set lower bound of all ireversible reactions to 0
    [m.set_constraint(k, lower_bound=0) for k, (lb, ub) in m.get_reactions_bounds().items() if lb > 0]

    # Close exachnge reactions bounds: lower = upper = 0
    [m.set_constraint(r, lower_bound=0, upper_bound=0) for r in m.get_exchanges()]

    # Run FVA
    fva = FVA(m)

    # Identify reactions in futile cycles
    futile_cycle_reactions = dict((r, 1) for r, (lb, ub) in fva.items() if lb < 0 or ub > 0)

    # Verbose
    if verbose > 0:
        print '[INFO] ' + str(len(futile_cycle_reactions)) + ' futile cycle reactions found: ', futile_cycle_reactions.keys()

    # Fix reactions values and minimise futile cycle reactions
    if len(futile_cycle_reactions) > 0:
        for i in samples.index:
            [m.set_constraint(r, samples.loc[i, r], samples.loc[i, r]) for r in samples.columns if r not in futile_cycle_reactions]

            fba = FBA(m, futile_cycle_reactions, maximize=False, opt_abs_values=True)

            for futile_r in futile_cycle_reactions:
                samples.ix[i, futile_r] = fba.values[futile_r]

    return samples

