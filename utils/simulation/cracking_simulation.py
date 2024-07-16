import random
import numpy as np
from matplotlib import pyplot as plt
import math




def crack_molecules(reactor,model):
    # determine number of zeros in each molecule
    no_zeros = [molecule.count(0) for molecule in reactor]
    if sum(no_zeros) == 0:
        return reactor
    #select a molecule with the probaility of the number of zeros
    selected_molecule = random.choices(reactor, weights=no_zeros)
    selected_list = [i for i, e in enumerate(selected_molecule[0]) if e == 0]
    #select a zero in the molecule based on the selected model
    normalized_probabilities = model(selected_molecule[0])

    selected_zero = random.choices(selected_list, weights=normalized_probabilities, k=1)
 
    #get the index of the selected molecule
    index = reactor.index(selected_molecule[0])
    #split the molecule at the selected zero, replace the orginal molecule with the two new molecules
    reactor[index] = selected_molecule[0][selected_zero[0]+1:]
    reactor.append(selected_molecule[0][:selected_zero[0]])
    return reactor

def remove_short_molecules(reactor, min_length):
    products = []
    reactor_remains = []
    for molecule in reactor:
        if sum(molecule) <= min_length:
            products.append(molecule)
        else:
            reactor_remains.append(molecule)
    return products, reactor_remains

def status(list):
    lenght = len(list)
    mass = sum([sum(molecule) for molecule in list])
    try: Mn = mass/lenght
    except ZeroDivisionError: Mn = 0
    return [lenght, mass, Mn]