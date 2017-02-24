try:
    import matplotlib.pyplot as plt
except:
    raise
import networkx as nx
G=nx.cycle_graph(24)
pos=nx.spring_layout(G,iterations=200)
print(len(G.edges()))
nx.draw(G,pos,node_color=[i for i in range(24)],node_size=800,cmap=plt.cm.jet, edge_color=[i for i in range(24)], edge_cmap=plt.cm.ocean)
plt.savefig("node_colormap.png") # save as png
plt.show() # display