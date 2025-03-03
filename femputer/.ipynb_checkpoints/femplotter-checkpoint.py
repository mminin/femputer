import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import pandas as pd

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

def plot_network_with_beams(nodes, elements, beam_assignments=None, displacement_multiplier=1.0, plot_deformed=False):
    """
    Plot the truss network with color-coded beams based on assigned shapes.

    Parameters:
        nodes (list): List of Node objects.
        elements (list): List of Element objects.
        beam_assignments (dict or None): Mapping of element IDs to beam shapes. If None, use material shapes from elements.
        displacement_multiplier (float): Factor for scaling displacements.
    """
    # Determine beam shapes and colors
    if beam_assignments:
        beam_shapes = list(set(beam_assignments.values()))
        color_map = {shape: cm.tab20(i / len(beam_shapes)) for i, shape in enumerate(beam_shapes)}
    else:
        # Use material definitions for shapes if no beam assignments
        beam_shapes = list(set(element.shape.key for element in elements))
        color_map = {shape: cm.tab20(i / len(beam_shapes)) for i, shape in enumerate(beam_shapes)}

    plt.figure(figsize=(12, 8))
    labels_shown = set()  # Track labels already added to the legend

    # Plot undeformed truss (blue lines)
    for element in elements:
        x_coords = [element.node1.x, element.node2.x]
        y_coords = [element.node1.y, element.node2.y]
        if 'Undeformed' not in labels_shown:
            plt.plot(x_coords, y_coords, 'b-', linewidth=1, label='Undeformed')
            labels_shown.add('Undeformed')
        else:
            plt.plot(x_coords, y_coords, 'b-', linewidth=1)

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

        # Determine beam shape
        if beam_assignments and element.id in beam_assignments:
            beam_shape = beam_assignments[element.id]
        else:
            beam_shape = element.shape.key  # Use material key as the shape

        # Get color for the beam shape
        color = color_map[beam_shape]

        if beam_shape not in labels_shown:
            plt.plot(x_coords, y_coords, color=color, linestyle='--', linewidth=2, label=beam_shape)
            labels_shown.add(beam_shape)
        else:
            plt.plot(x_coords, y_coords, color=color, linestyle='--', linewidth=2)

    # Plot nodes with proper label offset
    for node in nodes:
        plt.plot(node.x, node.y, 'ro')  # Undeformed position
        
        # Add labels with offset to prevent overlap
        offset = 0.5  # You can adjust the offset to be larger if needed
        #print(node.x,node.y)
        plt.text(node.x, node.y + offset+5, f" {node.id}", fontsize=10, verticalalignment='bottom', 
            horizontalalignment='right')

        # Deformed position
        if plot_deformed:
            deformed_x = node.x + displacement_multiplier * (node.displacement[0] or 0)
            deformed_y = node.y + displacement_multiplier * (node.displacement[1] or 0)
            plt.plot(deformed_x, deformed_y, 'go')

    # Add legend and labels
    plt.xlabel('X Coordinate')
    plt.ylabel('Y Coordinate')
    plt.title('Truss Network with Beam Assignments' if beam_assignments else 'Truss Network Using Material Shapes')
    
    # Adjust the axis to fit the entire grid
    plt.axis('equal')
    plt.grid(True)
    plt.tight_layout()  # Ensures the plot fits well in the figure area
    plt.legend(loc='upper left', fontsize='small')
    plt.show()


def plot_network_with_beams_20250117_1847(nodes, elements, beam_assignments=None, displacement_multiplier=1.0):
    """
    Plot the truss network with color-coded beams based on assigned shapes.

    Parameters:
        nodes (list): List of Node objects.
        elements (list): List of Element objects.
        beam_assignments (dict or None): Mapping of element IDs to beam shapes. If None, use material shapes from elements.
        displacement_multiplier (float): Factor for scaling displacements.
    """
    # Determine beam shapes and colors
    if beam_assignments:
        beam_shapes = list(set(beam_assignments.values()))
        color_map = {shape: cm.tab20(i / len(beam_shapes)) for i, shape in enumerate(beam_shapes)}
    else:
        # Use material definitions for shapes if no beam assignments
        beam_shapes = list(set(element.shape.key for element in elements))
        color_map = {shape: cm.tab20(i / len(beam_shapes)) for i, shape in enumerate(beam_shapes)}

    plt.figure(figsize=(12, 8))
    labels_shown = set()  # Track labels already added to the legend

    # Plot undeformed truss (blue lines)
    for element in elements:
        x_coords = [element.node1.x, element.node2.x]
        y_coords = [element.node1.y, element.node2.y]
        if 'Undeformed' not in labels_shown:
            plt.plot(x_coords, y_coords, 'b-', linewidth=1, label='Undeformed')
            labels_shown.add('Undeformed')
        else:
            plt.plot(x_coords, y_coords, 'b-', linewidth=1)

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

        # Determine beam shape
        if beam_assignments and element.id in beam_assignments:
            beam_shape = beam_assignments[element.id]
        else:
            beam_shape = element.shape.key  # Use material key as the shape

        # Get color for the beam shape
        color = color_map[beam_shape]

        if beam_shape not in labels_shown:
            plt.plot(x_coords, y_coords, color=color, linestyle='--', linewidth=2, label=beam_shape)
            labels_shown.add(beam_shape)
        else:
            plt.plot(x_coords, y_coords, color=color, linestyle='--', linewidth=2)

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
    plt.title('Truss Network with Beam Assignments' if beam_assignments else 'Truss Network Using Material Shapes')
    plt.axis('equal')
    plt.grid(True)
    plt.legend(loc='upper left', fontsize='small')
    plt.show()




def plot_principal_stress_with_tension_compression(nodes, elements, displacement_multiplier=1.0, show_forces = False):
    """
    Plot the truss network with color-coded elements based on principal stress magnitudes, using different
    line styles for tension and compression.

    Parameters:
        nodes (list): List of Node objects.
        elements (list): List of Element objects.
        displacement_multiplier (float): Factor for scaling displacements.
    """
    # Normalize stress values for color mapping
    if not show_forces:
        stresses = [abs(element.stress) for element in elements]
    else:
        stresses = [abs(element.internal_force) for element in elements]
    max_stress = max(stresses) if stresses else 1  # Avoid division by zero
    min_stress = min(stresses) if stresses else 0

    # Create a color map for stress values
    norm = plt.Normalize(vmin=min_stress, vmax=max_stress)
    cmap = cm.plasma

    fig, ax = plt.subplots(figsize=(12, 8))

    # Plot undeformed truss (gray lines)
    for element in elements:
        x_coords = [element.node1.x, element.node2.x]
        y_coords = [element.node1.y, element.node2.y]
        ax.plot(x_coords, y_coords, 'gray', linewidth=1, label='Undeformed' if element.id == 1 else "")

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

        # Map stress to a color
        if not show_forces:
            color = cmap(norm(abs(element.stress)))
            # Determine line style based on stress sign
            line_style = '-' if element.stress > 0 else '--'
        else:
            color = cmap(norm(abs(element.internal_force)))
            # Determine line style based on stress sign
            line_style = '-' if element.internal_force > 0 else '--'

        # Plot the deformed element with stress-based color and line style
        ax.plot(x_coords, y_coords, color=color, linestyle=line_style, linewidth=3)

    # Plot nodes
    for node in nodes:
        ax.plot(node.x, node.y, 'ro')  # Undeformed position
        ax.text(node.x, node.y, f" {node.id}", fontsize=10, verticalalignment='bottom', horizontalalignment='right')

        # Deformed position
        deformed_x = node.x + displacement_multiplier * (node.displacement[0] or 0)
        deformed_y = node.y + displacement_multiplier * (node.displacement[1] or 0)
        ax.plot(deformed_x, deformed_y, 'go')

    # Add color bar for stress
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.8)
    cbar.set_label('Principal Stress (Pa)', fontsize=12)

    # Add legend and labels
    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    if not show_forces:
        ax.set_title('Truss Network with Principal Stresses (Tension: Solid, Compression: Dashed)')
    else:
        ax.set_title('Truss Network with Internal Forces (Tension: Solid, Compression: Dashed)')
    ax.axis('equal')
    ax.grid(True)
    plt.show()


def plot_with_intrusion_20250115_2008(nodes, elements, displacement_multiplier=1.0, show_forces=False, intrusion_node_ids=None):
    """
    Plot the truss network with color-coded elements based on principal stress magnitudes, using different
    line styles for tension and compression.

    Parameters:
        nodes (list): List of Node objects.
        elements (list): List of Element objects.
        displacement_multiplier (float): Factor for scaling displacements.
        intrusion_node_ids (list or None): List of node IDs to highlight as intrusion nodes.
    """
    # Normalize stress values for color mapping
    if not show_forces:
        stresses = [abs(element.stress) for element in elements]
    else:
        stresses = [abs(element.internal_force) for element in elements]
    max_stress = max(stresses) if stresses else 1  # Avoid division by zero
    min_stress = min(stresses) if stresses else 0

    # Create a color map for stress values
    norm = plt.Normalize(vmin=min_stress, vmax=max_stress)
    cmap = cm.plasma

    fig, ax = plt.subplots(figsize=(12, 8))

    # Plot undeformed truss (gray lines)
    for element in elements:
        x_coords = [element.node1.x, element.node2.x]
        y_coords = [element.node1.y, element.node2.y]
        ax.plot(x_coords, y_coords, 'gray', linewidth=1, label='Undeformed' if element.id == 1 else "")

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

        # Map stress to a color
        if not show_forces:
            color = cmap(norm(abs(element.stress)))
            # Determine line style based on stress sign
            line_style = '-' if element.stress > 0 else '--'
        else:
            color = cmap(norm(abs(element.internal_force)))
            # Determine line style based on stress sign
            line_style = '-' if element.internal_force > 0 else '--'

        # Plot the deformed element with stress-based color and line style
        ax.plot(x_coords, y_coords, color=color, linestyle=line_style, linewidth=3)

    # Plot nodes
    for node in nodes:
        # Check if the node is in the list of intrusion nodes
        if node.id in intrusion_node_ids:
            ax.plot(node.x, node.y, 'ro')  # Intrusion node in red
        else:
            ax.plot(node.x, node.y, 'go')  # Regular node in green

        # Deformed position (if applicable)
        deformed_x = node.x + displacement_multiplier * (node.displacement[0] or 0)
        deformed_y = node.y + displacement_multiplier * (node.displacement[1] or 0)
        #ax.plot(deformed_x, deformed_y, 'ko')  # Deformed node in black

    # Add color bar for stress
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.8)
    cbar.set_label('Principal Stress (Pa)', fontsize=12)

    # Add legend and labels
    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    if not show_forces:
        ax.set_title('Truss Network with Principal Stresses (Tension: Solid, Compression: Dashed)')
    else:
        ax.set_title('Truss Network with Internal Forces (Tension: Solid, Compression: Dashed)')
    ax.axis('equal')
    ax.grid(True)
    plt.show()

def plot_with_intrusion(nodes, elements, displacement_multiplier=1.0, 
        show_forces=False, intrusion_node_ids=None, plot_deformed=False):
    """
    Plot the truss network with color-coded elements based on principal stress magnitudes, using different
    line styles for tension and compression.

    Parameters:
        nodes (list): List of Node objects.
        elements (list): List of Element objects.
        displacement_multiplier (float): Factor for scaling displacements.
        intrusion_node_ids (list or None): List of node IDs to highlight as intrusion nodes.
    """
    # Normalize stress values for color mapping
    if not show_forces:
        stresses = [abs(element.stress) for element in elements]
    else:
        stresses = [abs(element.internal_force) for element in elements]
    max_stress = max(stresses) if stresses else 1  # Avoid division by zero
    min_stress = min(stresses) if stresses else 0

    # Create a color map for stress values
    #norm = plt.Normalize(vmin=min_stress, vmax=max_stress)
    norm = plt.Normalize(vmin=min_stress, vmax=max_stress)
    cmap = cm.plasma

    fig, ax = plt.subplots(figsize=(12, 8))

    if plot_deformed:
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

            # Map stress to a color
            if not show_forces:
                color = cmap(norm(abs(element.stress)))
                # Determine line style based on stress sign
                line_style = '-' if element.stress > 0 else '--'
            else:
                color = cmap(norm(abs(element.internal_force)))
                # Determine line style based on stress sign
                line_style = '-' if element.internal_force > 0 else '--'

            # Plot the deformed element with stress-based color and line style
            ax.plot(x_coords, y_coords, color=color, linestyle=line_style, linewidth=3)
    else:
        for element in elements:
            x_coords = [
                element.node1.x or 0,
                element.node2.x or 0
            ]
            y_coords = [
                element.node1.y or 0,
                element.node2.y or 0
            ]

            # Map stress to a color
            if not show_forces:
                color = cmap(norm(abs(element.stress)))
                # Determine line style based on stress sign
                line_style = '-' if element.stress > 0 else '--'
            else:
                color = cmap(norm(abs(element.internal_force)))
                # Determine line style based on stress sign
                line_style = '-' if element.internal_force > 0 else '--'

            # Plot the deformed element with stress-based color and line style
            ax.plot(x_coords, y_coords, color=color, linestyle=line_style, linewidth=3)

    # Plot nodes
    print([node.id for node in nodes])
    print(intrusion_node_ids)
    for node in nodes:
        # Check if the node is in the list of intrusion nodes
        if node.id not in [str(_) for _ in intrusion_node_ids]:
            #ax.plot(node.x, node.y, 'go')  # Intrusion node in red
            pass
        else:
            ax.plot(node.x, node.y, 'ro')  # Regular node in green

        # Deformed position (if applicable)
        #deformed_x = node.x + displacement_multiplier * (node.displacement[0] or 0)
        #deformed_y = node.y + displacement_multiplier * (node.displacement[1] or 0)

    # Add color bar for stress
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.8)
    cbar.set_label('Principal Stress (Pa)', fontsize=12)

    # Add legend and labels
    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    if not show_forces:
        ax.set_title('Truss Network with Principal Stresses (Tension: Solid, Compression: Dashed)')
    else:
        ax.set_title('Truss Network with Internal Forces (Tension: Solid, Compression: Dashed)')
    ax.axis('equal')
    ax.grid(True)
    plt.show()




def printout_element_results(elements):
    for element in elements:
        print(f"Element {element.id} ({element.node1.id}-{element.node2.id}): "+
              f"{element.shape.key}, area {element.shape.area} m2, "+
              f"stress {round(element.stress/1000,2)} kPa, "+
              f"{round((element.stress/1000)/6895, 2)} ksi, "+
              f"forces {round((element.internal_force)/4.44822, 2)} lbs")


def printout_element_results_w_costs(elements, beam_assignments, pricing_table):
    """
    Print element results, including costs and stresses, and compute total cost.

    Parameters:
        elements (list): List of Element objects.
        beam_assignments (dict or None): Mapping of element IDs to beam shapes. If None, use material names.
        pricing_table (DataFrame): Pandas DataFrame containing beam prices and properties.
    """
    total_cost = 0
    for element in elements:
        # Determine the shape to use for pricing
        if beam_assignments and element.id in beam_assignments:
            beam_shape = beam_assignments[element.id]
        else:
            beam_shape = element.shape.key  # Use material name as the shape
        
        # Get the price per meter for the beam shape
        price_row = pricing_table.loc[pricing_table["Shape"] == beam_shape]
        if price_row.empty:
            raise ValueError(f"Price data for beam shape '{beam_shape}' not found in pricing table.")
        beam_price_per_meter = price_row["Cost per Meter (usd)"].values[0]
        
        # Calculate the length of the element (undeformed length)
        length = element.length
        
        # Compute the cost for this element
        element_cost = beam_price_per_meter * length
        total_cost += element_cost
        
        # Print element details
        print(f"Element {element.id} ({element.node1.id}-{element.node2.id}): " +
              f"{element.shape.key}, area {(element.shape.area*10000):.2f} cm2, " +
              f"stress {round(element.stress / 1000, 2)} kPa, " +
              f"{round((element.stress / 1000) / 6895, 2)} ksi, " +
              f"length {round(element.length, 2)} m, " +
              f"cost ${round(element_cost, 2)}, " +
              f"forces {round((element.internal_force)/4.44822, 2)} lbs")
    
    # Print total cost
    print(f"\nTotal Cost: ${round(total_cost, 2)}")

def make_df_element_results_w_costs(elements, beam_assignments, pricing_table):
    total_cost = 0
    data_list = []
    for element in elements:
        # Get the shape of the assigned beam
        beam_shape = beam_assignments[element.id]
        
        # Get the price per meter for the assigned beam
        beam_price_per_meter = pricing_table.loc[pricing_table["Shape"] == beam_shape, "Cost per Meter (usd)"].values[0]
        
        # Calculate the length of the element (undeformed length)
        length = element.length
        
        # Compute the cost for this element
        element_cost = beam_price_per_meter * length
        total_cost += element_cost
        element_descr = {"Element": {element.id}, "Nodes": f"{element.node1.id}-{element.node2.id}",
                    "material":element.shape.key, "area_cm2": element.shape.area*10000,
                    "internal_force_kN": element.internal_force / 1000,
                    "internal_force_lbs": (element.internal_force)/4.44822,
                    "stress_kPa": element.stress / 1000,
                    "stress_ksi": (element.stress / 1000) / 6895,
                    "length":element.length,
                    "cost": element_cost,
                    
                    }
        data_list.append(element_descr)
    return pd.DataFrame(data_list)




