import networkx as nx
from numpy.random import exponential
import time
N = 500

weights = exponential(1, N*N)

g = nx.complete_graph(N)

start = time.time()
for i in range(0, N):
    for j in range(0, N):
        g[i][j] = weights[i+j]

end = time.time()

print(end-start)
