#!/usr/bin/env python3
"""Random kbtree instance generator (binary CSPs whose constraint graph is a tree of cliques).

Parameters (w, s, h, d, t, k):
  w + 1  clique size
  s      separator size (variables shared by a clique and its parent)
  h      height of the clique tree (binary tree with h levels, 2^h - 1 cliques)
  d      domain size (values 0..d-1)
  t      tightness: percentage of the d*d tuples of each constraint that are forbidden (unit cost)
  k      number of instances

Every pair of variables of a clique gets one binary constraint whose floor(t*d*d/100) forbidden tuples are
drawn uniformly. Instances are written in XCSP 2.1 like the reference kbtree instances: relations listing
the forbidden tuples (semantics "conflicts"), or the allowed ones (semantics "supports") when they are fewer.

Output: <out>/kbtree-w-s-h-d/kbtree-w-s-h-d-t/normalized-kbtree-w-s-h-d-t-NN_ext.xml
For each instance, a variable ordering following the clique tree is written next to it (same name, .order):
the variables of the root cluster first, then a depth-first traversal where each cluster contributes its new
(non-separator) variables, so a whole subtree is assigned before the next one. One line per cluster.
"""

import argparse
import itertools
import os
import random


def generate(w, s, h, d, t, rng, branching=2):
    """Returns (nvars, constraints, order) where constraints is a sorted list of ((i, j), forbidden tuples) and
    order lists, for each cluster in depth-first order, the variables it introduces."""
    clique_size = w + 1
    root = list(range(clique_size))
    cliques = [root]
    new_vars = {0: root}   # cluster index -> variables it introduces
    children = {0: []}     # cluster index -> child cluster indices
    nvars = clique_size
    level = [0]
    for _ in range(h - 1):
        next_level = []
        for parent in level:
            for _ in range(branching):
                sep = rng.sample(cliques[parent], s)
                new = list(range(nvars, nvars + clique_size - s))
                nvars += clique_size - s
                child = len(cliques)
                cliques.append(sorted(sep + new))
                new_vars[child] = new
                children[child] = []
                children[parent].append(child)
                next_level.append(child)
        level = next_level

    edges = sorted({(i, j) for c in cliques for i, j in itertools.combinations(c, 2)})
    all_tuples = list(itertools.product(range(d), repeat=2))
    nconflicts = t * d * d // 100
    constraints = [(e, sorted(rng.sample(all_tuples, nconflicts))) for e in edges]

    # depth-first traversal of the clique tree, first child first
    order = []
    stack = [0]
    while stack:
        c = stack.pop()
        order.append(new_vars[c])
        stack.extend(reversed(children[c]))
    return nvars, constraints, order


def write_order(path, order):
    """Writes the variable ordering, one line per cluster (names separated by spaces)."""
    with open(path, "w") as f:
        f.write("# variable ordering following the clique tree (depth-first), one line per cluster\n")
        for vars_ in order:
            f.write(" ".join("V%d" % v for v in vars_) + "\n")


def write_xcsp(path, nvars, d, constraints):
    """Writes the instance in XCSP 2.1, sharing identical relations between constraints."""
    all_tuples = set(itertools.product(range(d), repeat=2))
    relations = {}  # (semantics, tuples) -> relation name
    refs = []
    for _, forbidden in constraints:
        if len(forbidden) > d * d / 2:
            key = ("supports", tuple(sorted(all_tuples - set(forbidden))))
        else:
            key = ("conflicts", tuple(forbidden))
        if key not in relations:
            relations[key] = "R%d" % len(relations)
        refs.append(relations[key])

    with open(path, "w") as f:
        f.write('<?xml version="1.0" encoding="UTF-8"?>\n<instance>\n')
        f.write('<presentation maxConstraintArity="2" format="XCSP 2.1"/>\n')
        f.write('<domains nbDomains="1">\n <domain name="D0" nbValues="%d">0..%d</domain>\n</domains>\n' % (d, d - 1))
        f.write('<variables nbVariables="%d">\n' % nvars)
        for i in range(nvars):
            f.write(' <variable name="V%d" domain="D0"/>\n' % i)
        f.write('</variables>\n<relations nbRelations="%d">\n' % len(relations))
        for (semantics, tuples), name in relations.items():
            f.write('<relation name="%s" arity="2" nbTuples="%d" semantics="%s">\n' % (name, len(tuples), semantics))
            f.write("|".join("%d %d" % tup for tup in tuples) + "\n</relation>\n")
        f.write('</relations>\n<constraints nbConstraints="%d">\n' % len(constraints))
        for c, ((i, j), _) in enumerate(constraints):
            f.write(' <constraint name="C%d" arity="2" scope="V%d V%d" reference="%s"/>\n' % (c, i, j, refs[c]))
        f.write('</constraints>\n</instance>\n')


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("w", type=int, help="clique size minus 1")
    p.add_argument("s", type=int, help="separator size")
    p.add_argument("h", type=int, help="height of the clique tree")
    p.add_argument("d", type=int, help="domain size")
    p.add_argument("t", type=int, help="tightness in percent (forbidden tuples)")
    p.add_argument("k", type=int, help="number of instances")
    p.add_argument("--out", default="Instance/kbtree", help="output root directory (default: Instance/kbtree)")
    p.add_argument("--seed", type=int, default=None, help="random seed (default: random)")
    p.add_argument("--branching", type=int, default=2, help="children per clique (default: 2)")
    a = p.parse_args()

    if not (0 <= a.s <= a.w):
        p.error("separator size s must satisfy 0 <= s <= w")
    if a.w < 1 or a.h < 1 or a.d < 1 or a.k < 1 or a.branching < 1:
        p.error("w, h, d, k and branching must be positive")
    if not (0 <= a.t <= 100):
        p.error("tightness t must be in 0..100")

    rng = random.Random(a.seed)
    family = "kbtree-%d-%d-%d-%d" % (a.w, a.s, a.h, a.d)
    name = "%s-%d" % (family, a.t)
    outdir = os.path.join(a.out, family, name)
    os.makedirs(outdir, exist_ok=True)
    width = max(2, len(str(a.k)))
    for n in range(1, a.k + 1):
        nvars, constraints, order = generate(a.w, a.s, a.h, a.d, a.t, rng, a.branching)
        path = os.path.join(outdir, "normalized-%s-%0*d_ext.xml" % (name, width, n))
        write_xcsp(path, nvars, a.d, constraints)
        write_order(path[:-len(".xml")] + ".order", order)
        print("%s: %d variables, %d constraints" % (path, nvars, len(constraints)))


if __name__ == "__main__":
    main()
