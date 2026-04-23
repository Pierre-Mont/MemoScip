#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import copy
import os
import sys
import itertools

assert len(sys.argv) == 3, "Please specify INPUT and OUTPUT filenames."

# does the WCSP have one or more constant terms
has_constant_term = False

# classe pour écriture à largeur de texte contrôlée
class WidthFile:
    maxcol = 80

    def __init__(self, *x, **k):
        # expecting filename, mode
        self.f = open(*x, **k)
        self.col = 0

    def write(self, x):
        lines = x.splitlines()
        if not lines:
            return
        # print "outputting", lines
        if (self.col + len(lines[0])) >= self.maxcol:
            self.f.write("\n")
            self.col = 0
        for ln in lines[:-1]:
            self.f.write(ln + '\n')
        self.f.write(lines[-1])
        if len(lines) > 1:
            self.col = len(lines[-1])
        else:
            self.col += len(lines[-1])

    def flush(self):
        self.f.flush()

    def close(self):
        self.f.close()

    def __getattr__(self, name):
        return getattr(self.f, name)


# le nom des variables pour l'encodage des valeurs des domaines. les
# variables booléennes restent booléennes mais il faut compter les
# littéraux négatifs par ailleurs (format lp ne gère pas (1-x))
def domain_var(n, v):
    return " d%i_%i " % (n, v)

def mdomain_var(coeff, n, v):
    return "%+i d%i_%i " % (coeff, n, v)


# le nom des variables pour l'encodage des autres tuples
def tuple_var(tvar, tval):
    tvarval = list(zip(tvar, tval))
    # normalize tuple
    st = sorted(tvarval, key=lambda x: x[0])
    name = "t"
    for x in st:
        name = name + ("_%i_%i" % x)
    return name

# le produit cartésien des séquences (stockées dans une séquence vlist).
def product(vlist):
    return itertools.product(*vlist)

# enumerate all "tuples" on tvar (for var, if it appears in tvar, a
# single value val is used instead of thh full domain) generating the
# set of support tuples.
def enum_tuples(tvar, var, val):
    return product([ [val] if (var == ovar) else range(domains[ovar]) for ovar in tvar ])

# reading numbers
def read_num_vec(toks):
    return [int(x) for x in toks]

def read_int_tok(tok_iter):
    return int(tok_iter(1)[0])

# lire une définition de cost function. The cost table is a tuple based dictionary
def read_fun(tok_iter):
    n_var = read_int_tok(tok_iter)
    vars_ = read_num_vec(tok_iter(n_var))
    defcost = read_int_tok(tok_iter)
    if defcost == -1:
        defcost = tok_iter(1)[0]
        n_spec = 1
    else:
        n_spec = read_int_tok(tok_iter)
    tvo = sorted(list(zip(vars_, range(len(vars_)))), key=lambda x: x[0])
    ovars = tuple(x[0] for x in tvo)
    varorder = tuple(x[1] for x in tvo)
    specs = dict()
    for i in range(n_spec):
        if defcost != 'knapsackp':
            tc = read_num_vec(tok_iter(n_var + 1))
            if isinstance(defcost, str):
                specs = tc
            else:
                specs[tuple(tc[i] for i in varorder)] = tc[-1]
        else:
            Weight = []
            Weight.append(read_int_tok(tok_iter))
            for j in range(n_var):
                nbval = read_int_tok(tok_iter)
                Weight.append(nbval)
                for k in range(nbval):
                    Weight.append(read_int_tok(tok_iter))
                    Weight.append(read_int_tok(tok_iter))
            specs = Weight
    if isinstance(defcost, str):
        return vars_, defcost, specs
    else:
        return ovars, defcost, specs

# parcourir une cost function table
def iter_fun(vars_, defcost, specs):
    vardom = [range(domains[v]) for v in vars_]
    for t in itertools.product(*vardom):
        if t in specs:
            yield t, specs[t]
        else:
            yield t, defcost

# parcourir une cost function table en évitant les tuples d'un coût
# donné si possible (defcost)
def iter_funavoid(vars_, defcost, specs, avoid):
    if (defcost == avoid):
        for t in specs:
            yield t, specs[t]
    else:
        vardom = [range(domains[v]) for v in vars_]
        for t in itertools.product(*vardom):
            if t in specs:
                yield t, specs[t]
            else:
                yield t, defcost

# -------------  MAIN ---------------------

def token_iter(filename):
    with open(filename, 'r') as fh:
        for l in fh:
            for stok in l.strip().split(" "):
                for ttok in stok.strip().split("\t"):
                    if ttok:
                        yield ttok

tokens = token_iter(sys.argv[1])


def next_tokens(n):
    return [next(tokens) for i in range(n)]

# line_iter = open(sys.argv[1]).xreadlines()
output = WidthFile(sys.argv[2], 'w')

print("File %s opened" % sys.argv[1])

# reading parameters
name = next(tokens)
n_var, max_domain_size, n_fun, upper_bound = read_num_vec(next_tokens(4))

domains = read_num_vec(next_tokens(n_var))
n_fun = int(n_fun)
ub = int(upper_bound)
output.write("Minimize \n")
all_fun = [read_fun(next_tokens) for i in range(n_fun)]

# print "\nCost functions read."

# Output the criteria. Do not integrate zero or "infinite" cost
# components here. Zero is useless, "infinite" will be handled as
# linear constraints
negative_litterals = 0

for vars_, defcost, specs in all_fun:
    if isinstance(defcost, str):
        continue
    n_vars = len(vars_)
    if (n_vars == 0):
        has_constant_term = 1
        output.write(' +%i t ' % defcost)
    else:
        for t, cost in iter_funavoid(vars_, defcost, specs, 0):
            if cost == 0:
                continue
            if n_vars == 1:
                output.write(mdomain_var(cost, vars_[0], t[0]))
                if (t[0] == 1 and domains[vars_[0]] <= 2):
                    negative_litterals = negative_litterals + cost
            else:
                output.write(' +%i %s ' % (cost, tuple_var(vars_, t)))
if negative_litterals and 1 == 0:
    has_constant_term = 1
    output.write(" +%i t" % negative_litterals)

# print "Criteria generated."

output.write("\n\nSubject to\n\n")

# Set of tuple vars that need not be used

ub_tuplevars = set()

# Hard constraints: for every value with cost >=ub, we forbid it
# explicitely. Tuples variables are just removed.
nbslack = 0
nblin=0
nbdom=0
for vars_, defcost, specs in all_fun:
    if isinstance(defcost, str):
        addtocap = 0
        if defcost == 'knapsack':
            output.write("linear_constraint_%i: " % nblin)
            for i, v in enumerate(vars_):
                if int(specs[i+1] >= 0):
                    output.write('+%i d%i_1' % (specs[i+1], v))
                else:
                    output.write('+%i d%i_0' % (-specs[i+1], v))
                    addtocap = addtocap - specs[i+1]
            output.write(' >= %i\n\n' % (specs[0] + addtocap,))
            nblin = nblin + 1
        elif defcost == 'knapsackp':
            last = 1
            output.write("linear_constraint_%i: " % nblin)
            for i, v in enumerate(vars_):
                nbval = specs[last]
                if (domains[v] == 1):
                    val1 = 0
                    val0 = 0
                    for j in range(nbval):
                        if int(specs[last+1+2*j]) == 1:
                            val1 = int(specs[last+2+2*j])
                        else:
                            val0 = int(specs[last+2+2*j])
                    if val1 >= 0:
                        output.write('+%i d%i_1' % (val1, v))
                    else:
                        output.write('%i d%i_1' % (val1, v))
                    if val0 >= 0:
                        output.write('+%i d%i_0' % (val0, v))
                    else:
                        output.write('%i d%i_0' % (val0, v))
                else:
                    for j in range(nbval):
                        coeff = specs[last+2+2*j]
                        idx = specs[last+1+2*j]
                        if int(coeff >= 0):
                            output.write('+%i d%i_%i' % (coeff, v, idx))
                        else:
                            output.write('%i d%i_%i' % (coeff, v, idx))
                last = last + nbval*2 + 1
            output.write(' >= %i\n\n' % (specs[0],))
            nblin = nblin + 1
        continue
    n_vars = len(vars_)
    for t, cost in iter_funavoid(vars_, defcost, specs, 0):
        if cost < ub:
            continue
        if n_vars == 1:
            output.write('%s = %i\n\n' % (mdomain_var(1, vars_[0], t[0]), -(domains[vars_[0]] == 2 and t[0] == 1)))

# print "Hard constraints generated."

# Direct encoding. Exactly one value constraint. Boolean variables are not included here.
for i, dom in enumerate(domains):
    if (dom > 1):
        output.write("domain_constraint_%i: " % nbdom)
        for v in range(dom):
            output.write(mdomain_var(1, i, v))
        output.write(" = 1\n\n")
        nbdom = nbdom + 1
    if (dom == 1):
        output.write("domain_constraint_%i: " % nbdom)
        output.write("%s = 1\n\n" % domain_var(i, 0))
        nbdom = nbdom + 1

# print "Domain constraints generated."

# marginal consistency: one value selected iff one associated tuple selected.
# if several functions have the same cost, we need only to do it once
scopes = set(vars_ for vars_, defcost, specs in (f for f in all_fun if len(f[0]) >= 2 and not isinstance(f[1], str)))
# print "%i different scopes detected." % len(scopes)

for vars_ in scopes:
    for va in vars_:
        for a in range(domains[va]):
            for b in enum_tuples(vars_, va, a):
                tname = tuple_var(vars_, b)
                if tname not in ub_tuplevars:
                    output.write("+1 %s " % tname)
            output.write(" %s " % mdomain_var(-1, va, a))
            output.write("= 0\n\n")

# print "Marginal consistency constraints generated."

if has_constant_term:
    output.write("t = 1\n\n")

# print "Tuple bounds generated."

output.write("\n\nBinary\n\n")
# indicate 0/1 variables (direct encoding).
for i, dom in enumerate(domains):
    if (dom > 1):
        for v in range(dom):
            output.write("%s " % domain_var(i, v))
    else:
        output.write("%s " % domain_var(i, 0))
for vars_, defcost, specs in (f for f in all_fun if len(f[0]) >= 2):
    if isinstance(defcost, str):
        continue
    for b in enum_tuples(vars_, -1, -1):
        tname = tuple_var(vars_, b)
        if tname not in ub_tuplevars:
            output.write("%s " % tname)

if has_constant_term:
    output.write("t")

output.write("\n\nEnd")
# print "Domain binaries generated."
print("File " + sys.argv[2] + " has been created")