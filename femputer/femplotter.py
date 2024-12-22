import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

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
    plt.title(f'FEM Network Visualization: Undeformed and Deformed Shapes, x{displacement_multiplier}')
    plt.axis('equal')
    plt.grid(True)
    #plt.legend()
    plt.show()

def plot_network_with_beams(nodes, elements, beam_assignments=None, displacement_multiplier=1.0):
    """
    Plot the truss network with color-coded beams based on assigned shapes.

    Parameters:
        nodes (list): List of Node objects.
        elements (list): List of Element objects.
        beam_assignments (dict or None): Mapping of element IDs to beam shapes. If None, no beam colors are applied.
        displacement_multiplier (float): Factor for scaling displacements.
    """
    # Assign colors based on beam shapes or default to one color if no assignments
    if beam_assignments:
        beam_shapes = list(set(beam_assignments.values()))
        color_map = {shape: cm.tab20(i / len(beam_shapes)) for i, shape in enumerate(beam_shapes)}
    else:
        beam_shapes = []
        color_map = None

    plt.figure(figsize=(12, 8))

    # Plot undeformed truss (blue lines)
    for element in elements:
        x_coords = [element.node1.x, element.node2.x]
        y_coords = [element.node1.y, element.node2.y]
        plt.plot(x_coords, y_coords, 'b-', linewidth=1, label='Undeformed' if element.id == 1 else None)

    # Plot deformed truss with color-coded elements
    for element in elements:
        x_coords = [
            element.node1.x + displacement_multiplier * (element.node1.displacement[0] or 0),
            element.node2.x + displacement_multiplier * (element.node2.displacement[0] or 0),
        ]
        y_coords = [
            element.node1.y + displacement_multiplier * (element.node1.displacement[1] or 0),
            element.node2.y + displacement_multiplier * (element.node2.displacement[1] or 0),
        ]

        if beam_assignments:
            beam_shape = beam_assignments[element.id]
            color = color_map[beam_shape]
            plt.plot(x_coords, y_coords, color=color, linestyle='--', linewidth=2, label=beam_shape if element.id == 1 else "")
        else:
            plt.plot(x_coords, y_coords, 'r--', linewidth=2)  # Default color for all beams

    # Plot nodes
    for node in nodes:
        plt.plot(node.x, node.y, 'ro')  # Undeformed position
        plt.text(node.x, node.y, f" {node.id}", fontsize=10, verticalalignment='bottom', horizontalalignment='right')

        # Deformed position
        deformed_x = node.x + displacement_multiplier * (node.displacement[0] or 0)
        deformed_y = node.y + displacement_multiplier * (node.displacement[1] or 0)
        plt.plot(deformed_x, deformed_y, 'go')

    # Add legend and labels
    plt.xlabel('X Coordinate')
    plt.ylabel('Y Coordinate')
    plt.title('Truss Network with Beam Assignments' if beam_assignments else 'Truss Network Without Beam Assignments')
    plt.axis('equal')
    plt.grid(True)
    if beam_assignments:
        plt.legend(loc='upper left', fontsize='small')
    plt.show()
















