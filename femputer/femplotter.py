import matplotlib.pyplot as plt

# Plot the network using Matplotlib
def plot_network(nodes, elements, displacement_multiplier=1.0):
    plt.figure(figsize=(10, 6))
    
    # Plot undeformed network
    for element in elements:
        x_coords = [element.node1.x, element.node2.x]
        y_coords = [element.node1.y, element.node2.y]
        plt.plot(x_coords, y_coords, 'b-', linewidth=2, label='Undeformed' if element.id == 0 else "")
    
    # Plot deformed network
    for element in elements:
        x_coords = [
            element.node1.x + displacement_multiplier * element.node1.displacement[0],
            element.node2.x + displacement_multiplier * element.node2.displacement[0]
        ]
        y_coords = [
            element.node1.y + displacement_multiplier * element.node1.displacement[1],
            element.node2.y + displacement_multiplier * element.node2.displacement[1]
        ]
        plt.plot(x_coords, y_coords, 'r--', linewidth=2, label='Deformed' if element.id == 0 else "")
    
    # Plot nodes (undeformed)
    for node in nodes:
        plt.plot(node.x, node.y, 'ro')
        plt.text(node.x, node.y, f" {node.id}", fontsize=12, verticalalignment='bottom', horizontalalignment='right')
    
    plt.xlabel('X Coordinate')
    plt.ylabel('Y Coordinate')
    plt.title('FEM Network Visualization: Undeformed and Deformed Shapes')
    plt.axis('equal')
    plt.grid(True)
    #plt.legend()
    plt.show()
