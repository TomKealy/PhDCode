import networkx as nx
filename = "net_data.dat"
P = 50
d = 0.25
G = nx.random_geometric_graph(P,d)
Adj = nx.adj_matrix(G)
pos = nx.get_node_attributes(G,'pos')
if nx.is_connected(G):
	conn = 1
else:
	conn = 0

FILE = open(filename, "w")
FILE.write('conn = '+str(int(conn))+'\n')
FILE.write('Adj = \n')
for i_1 in range(0, P):
	for i_2 in range(0, P):
		num = int(Adj[i_1,i_2])
		FILE.write(str(num) + ' ')
	FILE.write('\n')
FILE.write('pos = \n')
for i in range(0, P):
    FILE.write(str(pos[i]) + '\n')
FILE.write('\n')
FILE.close()
