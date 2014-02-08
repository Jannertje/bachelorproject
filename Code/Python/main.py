from algos import Greedy as Algo
from sorting import Sort as s

f = lambda x: 10
inst = Algo(f, s)
inst.iteration(leaves = [[0,0.25],[0.25,0.5],[0.5,1.0]])
