#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy
from concurrent.futures import ProcessPoolExecutor
from collections import namedtuple, Counter
from functools import partial

Column = namedtuple('Column', ('name', 'data'))
Result = namedtuple('Result', ('name', 'classes', 'entropy'))

def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('--calls', required=True)

    parser.add_argument('--replicates', default=100, type=int)

    parser.add_argument('--cores', default=1, type=int)

    return parser.parse_args()

def observations(calls):

    for name in calls:
        #print(name)
        data = [x for x in calls[name] if x > 0]
        yield Column(name, data)

def shannon(s):

    assert(type(s) in (str, list, tuple))

    length = len(s)

    classes = Counter(s)

    def f(x, log=numpy.log2):
        a = classes[x]
        return (a * log(a / length)) / length

    return -sum(f(x) for x in classes)

def simulate_discovery(name, col, replicates):

    discrete_class_results, shannon_results = [], []
    mean = numpy.mean  # optimize dereference

    for rep in range(replicates):
        numpy.random.seed(rep)  # reproducibility

        cur = [x for x in col]
        numpy.random.shuffle(cur)

        discrete_classes = []
        shannons = []

        for i, element in enumerate(cur, 1):

            observed = cur[:i]
            discrete_classes.append( len(set(observed)) )

            #shannons.append( shannon(observed) )

        discrete_class_results.append(discrete_classes)
        #shannon_results.append(shannons)

    return name, [mean(x) for x in zip(*discrete_class_results)], None
             #[mean(y) for y in zip(*shannon_results)]

def simulation_worker(name, data,  replicates):
    name, classes, entropy = simulate_discovery(name, data, replicates)
    return Result(name, classes, entropy)

def discovery_probabilities(calls, replicates, cores):

    worker = partial(simulation_worker, replicates=replicates)

    with ProcessPoolExecutor(cores) as ppe:

        results = ppe.map(worker, *zip(*observations(calls)))

    return sorted(results, key=lambda x: x.name)

def main():

    args = arguments()

    calls = pd.read_csv(args.calls, index_col='genomes')

    discovery_probabilities(calls, args.replicates, args.cores)

if __name__ == '__main__':
    main()
