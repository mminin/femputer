import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd

from mpl_toolkits.mplot3d import Axes3D
import ipywidgets as widgets
from IPython.display import display


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
        offset = 0  # You can adjust the offset to be larger if needed
        #print(node.x,node.y)
        plt.text(node.x, node.y + offset, f" {node.id}", fontsize=10, verticalalignment='bottom', 
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



def plot_network_with_beams_3d(nodes, elements, beam_assignments=None, displacement_multiplier=1.0, plot_deformed=False):
    """
    Interactive 3D plot of the truss network with color-coded beams based on assigned shapes.

    Parameters:
        nodes (list): List of Node objects.
        elements (list): List of Element objects.
        beam_assignments (dict or None): Mapping of element IDs to beam shapes. If None, use material shapes from elements.
        displacement_multiplier (float): Factor for scaling displacements.
        plot_deformed (bool): Whether to plot the deformed shape of the truss.
    """
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.view_init(elev=30, azim=30)  # Set an initial viewing angle
    
    # Create interactive sliders for elevation and azimuth
    elev_slider = widgets.FloatSlider(min=0, max=180, step=1, value=30, description='Elevation')
    azim_slider = widgets.FloatSlider(min=0, max=360, step=1, value=30, description='Azimuth')
    
    def update_view(elev, azim):
        ax.view_init(elev=elev, azim=azim)
        plt.draw()
    
    widgets.interactive(update_view, elev=elev_slider, azim=azim_slider)
    display(elev_slider, azim_slider)

    # Determine beam shapes and colors
    if beam_assignments:
        beam_shapes = list(set(beam_assignments.values()))
        color_map = {shape: cm.tab20(i / len(beam_shapes)) for i, shape in enumerate(beam_shapes)}
    else:
        beam_shapes = list(set(element.shape.key for element in elements))
        color_map = {shape: cm.tab20(i / len(beam_shapes)) for i, shape in enumerate(beam_shapes)}
    
    labels_shown = set()

    # Plot undeformed truss (blue lines)
    for element in elements:
        x_coords = [element.node1.x, element.node2.x]
        y_coords = [element.node1.y, element.node2.y]
        z_coords = [element.node1.z, element.node2.z]
        
        if 'Undeformed' not in labels_shown:
            ax.plot(x_coords, y_coords, z_coords, 'b-', linewidth=1, label='Undeformed')
            labels_shown.add('Undeformed')
        else:
            ax.plot(x_coords, y_coords, z_coords, 'b-', linewidth=1)

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
        z_coords = [
            element.node1.z + displacement_multiplier * (element.node1.displacement[2] or 0),
            element.node2.z + displacement_multiplier * (element.node2.displacement[2] or 0),
        ]
        
        beam_shape = beam_assignments[element.id] if beam_assignments and element.id in beam_assignments else element.shape.key
        color = color_map[beam_shape]

        if beam_shape not in labels_shown:
            ax.plot(x_coords, y_coords, z_coords, color=color, linestyle='--', linewidth=2, label=beam_shape)
            labels_shown.add(beam_shape)
        else:
            ax.plot(x_coords, y_coords, z_coords, color=color, linestyle='--', linewidth=2)

    # Plot nodes with proper label offset
    for node in nodes:
        ax.scatter(node.x, node.y, node.z, color='red', marker='o')
        ax.text(node.x, node.y, node.z, f" {node.id}", fontsize=10, verticalalignment='bottom', horizontalalignment='right')

    # Deformed position
    if plot_deformed:
        for node in nodes:
            deformed_x = node.x + displacement_multiplier * (node.displacement[0] or 0)
            deformed_y = node.y + displacement_multiplier * (node.displacement[1] or 0)
            deformed_z = node.z + displacement_multiplier * (node.displacement[2] or 0)
            ax.scatter(deformed_x, deformed_y, deformed_z, color='green', marker='o')
    
    # Labels and title
    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    ax.set_zlabel('Z Coordinate')
    ax.set_title('Interactive 3D Truss Network with Beam Assignments' if beam_assignments else 'Interactive 3D Truss Network')
    ax.legend(loc='upper left', fontsize='small')
    plt.show()


def plot_principal_stress_with_tension_compression_3d(nodes, elements, displacement_multiplier=1.0, show_forces=False, show_nodes=False, show_node_ids=False):
    """
    Interactive 3D plot of the truss network with principal stress visualization.
    Uses different colors and line styles for tension and compression elements.

    Parameters:
        nodes (list): List of Node objects.
        elements (list): List of Element objects.
        displacement_multiplier (float): Factor for scaling displacements.
        show_forces (bool): Whether to visualize internal forces instead of stress.
        show_nodes (bool): Whether to plot nodes.
        show_node_ids (bool): Whether to label node IDs.
    """
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.view_init(elev=30, azim=30)  # Initial view angle
    
    # Extract stress or force values for visualization
    values = [element.internal_force if show_forces else element.stress for element in elements]
    
    positive_stresses = [val for val in values if val > 0]
    negative_stresses = [val for val in values if val < 0]

    if positive_stresses:
        norm_pos = mcolors.Normalize(vmin=min(positive_stresses), vmax=max(positive_stresses))
    else:
        norm_pos = None
    
    if negative_stresses:
        norm_neg = mcolors.Normalize(vmin=min(negative_stresses), vmax=max(negative_stresses))
    else:
        norm_neg = None
    
    cmap_pos = cm.YlOrRd  # Yellow to Red for tension
    cmap_neg = cm.winter  # Blue shades for compression
    
    # Plot elements with color-coded stress
    for element in elements:
        x_coords = [
            element.node1.x + displacement_multiplier * (element.node1.displacement[0] or 0),
            element.node2.x + displacement_multiplier * (element.node2.displacement[0] or 0)
        ]
        y_coords = [
            element.node1.y + displacement_multiplier * (element.node1.displacement[1] or 0),
            element.node2.y + displacement_multiplier * (element.node2.displacement[1] or 0)
        ]
        z_coords = [
            element.node1.z + displacement_multiplier * (element.node1.displacement[2] or 0),
            element.node2.z + displacement_multiplier * (element.node2.displacement[2] or 0)
        ]
        
        value = element.internal_force if show_forces else element.stress
        
        if value == 0:
            color = 'black'  # Zero stress elements
            line_style = '-'
        elif value > 0:
            color = cmap_pos(norm_pos(value)) if norm_pos else 'red'
            line_style = '-'
        else:
            color = cmap_neg(norm_neg(value)) if norm_neg else 'blue'
            line_style = '--'
        
        ax.plot(x_coords, y_coords, z_coords, color=color, linestyle=line_style, linewidth=2)
    
    # Plot nodes if requested
    if show_nodes:
        for node in nodes:
            ax.scatter(node.x, node.y, node.z, color='red', marker='o')
            if show_node_ids:
                ax.text(node.x, node.y, node.z, f" {node.id}", fontsize=10, verticalalignment='bottom', horizontalalignment='right')
    
    # Labels and title
    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    ax.set_zlabel('Z Coordinate')
    ax.set_title('3D Principal Stress Distribution' if not show_forces else '3D Internal Force Distribution')
    
    plt.colorbar(cm.ScalarMappable(norm=norm_pos, cmap=cmap_pos), ax=ax, label='Tension')
    plt.colorbar(cm.ScalarMappable(norm=norm_neg, cmap=cmap_neg), ax=ax, label='Compression')
    plt.show()



def plot_principal_stress_with_tension_compression(nodes, elements, displacement_multiplier=1.0, show_forces=False, show_nodes=False, show_node_ids=False):
    """
    Plot the truss network with color-coded elements based on principal stress magnitudes, using different
    line styles for tension and compression, and black for elements with zero stress.

    Parameters:
        nodes (list): List of Node objects.
        elements (list): List of Element objects.
        displacement_multiplier (float): Factor for scaling displacements.
    """
    # Normalize stress values for color mapping
    if not show_forces:
        stresses = [element.stress for element in elements]
    else:
        stresses = [element.internal_force for element in elements]
    
    #print(max_stress, min_stress)
    positive_stresses = [stress for stress in stresses if stress>0]
    negative_stresses = [stress for stress in stresses if stress<0]

    #positive_max_stress = max(positive_stresses) if positive_stresses else 0  # Avoid division by zero
    #positive_min_stress = min(positive_stresses) if positive_stresses else 0

    #negative_max_stress = max(negative_stresses) if negative_stresses else 0  # Avoid division by zero
    #negative_min_stress = min(negative_stresses) if negative_stresses else 0

    # Ensure normalization range is valid
    if positive_stresses:
        positive_min_stress = min(positive_stresses)
        positive_max_stress = max(positive_stresses)
        if positive_min_stress == positive_max_stress:
            positive_min_stress -= 1
            positive_max_stress += 1
        norm_pos = mcolors.Normalize(vmin=positive_min_stress, vmax=positive_max_stress)
    else:
        norm_pos = None  # No positive stresses

    if negative_stresses:
        negative_min_stress = min(negative_stresses)
        negative_max_stress = max(negative_stresses)
        if negative_min_stress == negative_max_stress:
            negative_min_stress -= 1
            negative_max_stress += 1
        norm_neg = mcolors.Normalize(vmin=negative_min_stress, vmax=negative_max_stress)
    else:
        norm_neg = None  # No negative stresses

    #print("POSITIVE STRESSES", positive_max_stress, positive_min_stress)
    #print("NEGATIVE STRESSES", negative_max_stress, negative_min_stress)

    # Create two separate color maps: one for positive and one for negative values
    norm_pos = mcolors.Normalize(vmin=positive_min_stress, vmax=positive_max_stress)
    norm_neg = mcolors.Normalize(vmin=negative_min_stress, vmax=negative_max_stress)



    cmap_pos = cm.YlOrRd  # Yellow to Red for positive stresses (tension)
    cmap_neg = cm.winter  # Green to Blue for negative stresses (compression)

    fig, ax = plt.subplots(figsize=(12, 8))

    # Plot undeformed truss (gray lines)
    for element in elements:
        x_coords = [element.node1.x, element.node2.x]
        y_coords = [element.node1.y, element.node2.y]
        ax.plot(x_coords, y_coords, 'gray', linewidth=0.5, label='Undeformed' if element.id == 1 else "")

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

        if not show_forces:
            # Check for zero stress and plot black for zero stress elements
            if element.stress == 0:
                color = 'black'
                line_style = '-'  # Solid line for zero stress
            elif element.stress > 0:
                color_norm_pos = norm_pos(element.stress)
                print("color_norm_pos", color_norm_pos)
                color = cmap_pos(color_norm_pos)  # Yellow to Red for tension
                line_style = '-'  # Solid for tension
            else:
                color = cmap_neg(norm_neg(element.stress))  # Green to Blue for compression
                line_style = '-'  # Dashed for compression
        else:
            # Check for zero stress and plot black for zero stress elements
            if element.internal_force == 0:
                color = 'black'
                line_style = '-'  # Solid line for zero stress
            elif element.internal_force > 0:
                color_norm_pos = norm_pos(element.internal_force)
                print("color_norm_pos", color_norm_pos)
                color = cmap_pos(color_norm_pos)  # Yellow to Red for tension
                line_style = '-'  # Solid for tension
            else:
                color = cmap_neg(norm_neg(element.internal_force))  # Green to Blue for compression
                line_style = '-'  # Dashed for compression


        # Plot the deformed element with stress-based color and line style
        #print(element.stress, color)
        ax.plot(x_coords, y_coords, color=color, linestyle=line_style, linewidth=2)

    # Plot nodes
    if show_nodes:
        for node in nodes:
            ax.plot(node.x, node.y, 'ro')  # Undeformed position
            if show_node_ids:
                ax.text(node.x, node.y, f" {node.id}", fontsize=10, verticalalignment='bottom', horizontalalignment='right')

            # Deformed position
            deformed_x = node.x + displacement_multiplier * (node.displacement[0] or 0)
            deformed_y = node.y + displacement_multiplier * (node.displacement[1] or 0)
            ax.plot(deformed_x, deformed_y, 'go')

    # Add color bars for tension and compression
    sm_pos = cm.ScalarMappable(cmap=cmap_pos, norm=norm_pos)
    sm_pos.set_array([])
    cbar_pos = fig.colorbar(sm_pos, ax=ax, shrink=0.8)
    cbar_pos.set_label('Tension', fontsize=12)

    sm_neg = cm.ScalarMappable(cmap=cmap_neg, norm=norm_neg)
    sm_neg.set_array([])
    cbar_neg = fig.colorbar(sm_neg, ax=ax, shrink=0.8)
    cbar_neg.set_label('Compression', fontsize=12)

    # Add legend and labels
    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    if not show_forces:
        ax.set_title('Truss Network with Principal Stresses')
    else:
        ax.set_title('Truss Network with Internal Forces')
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
    #print([node.id for node in nodes])
    #print(intrusion_node_ids)
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


def plot_forces_on_nodes(nodes, forces, scale_factor=1e-6):
    """
    Plot the forces applied to the nodes as arrows.
    
    Parameters:
        nodes (list): List of Node objects.
        forces (dict): Dictionary of forces applied to each node.
        scale_factor (float): Factor to scale the arrows for better visibility.
    """
    # Extract the x and y coordinates of the nodes
    x_coords = [node.x for node in nodes]
    y_coords = [node.y for node in nodes]
    
    # Extract the force components for each node
    fx = [forces[str(i)]["x"] for i in range(len(nodes))]
    fy = [forces[str(i)]["y"] for i in range(len(nodes))]

    # Create the plot
    plt.figure(figsize=(12, 8))
    
    # Plot the nodes as dots
    plt.scatter(x_coords, y_coords, color='green', s=10, label='Nodes')

    # Plot the forces as arrows
    plt.quiver(x_coords, y_coords, fx, fy, angles='xy', scale_units='xy', scale=scale_factor, color='blue', label='Forces')

    # Set labels and title
    plt.xlabel('X Coordinate')
    plt.ylabel('Y Coordinate')
    plt.title('Force Distribution on Nodes')
    plt.grid(True)
    plt.legend()
    
    # Show the plot
    plt.axis('equal')
    plt.show()


def calculate_reaction_forces(nodes, elements, boundary_conditions):
    """
    Compute the total reaction forces at fixed nodes.
    
    Parameters:
        nodes (list): List of Node objects.
        elements (list): List of Element objects.
        boundary_conditions (dict): Dictionary defining fixed degrees of freedom for each node.

    Returns:
        dict: Total reaction force at each fixed node.
    """
    reaction_forces = {node.id: np.array([0.0, 0.0, 0.0]) for node in nodes if node.id in boundary_conditions}
    
    for element in elements:
        force = element.internal_force
        direction = np.array([
            (element.node2.x - element.node1.x),
            (element.node2.y - element.node1.y),
            (element.node2.z - element.node1.z)
        ])
        length = np.linalg.norm(direction)
        direction /= length  # Normalize direction vector
        
        force_vector = force * direction
        
        if element.node1.id in reaction_forces:
            reaction_forces[element.node1.id] += force_vector
        if element.node2.id in reaction_forces:
            reaction_forces[element.node2.id] -= force_vector
    
    return reaction_forces

import numpy as np

def get_boundary_conditions(nodes):
    """
    Extracts boundary conditions from node constraints.
    
    Parameters:
        nodes (list): List of Node objects.
    
    Returns:
        dict: Dictionary of fixed nodes and their constrained DOFs.
    """
    return {node.id: node.fixed for node in nodes if any(node.fixed.values())}

def calculate_reaction_forces(nodes, elements):
    """
    Compute the total reaction forces at fixed nodes.
    
    Parameters:
        nodes (list): List of Node objects.
        elements (list): List of Element objects.

    Returns:
        dict: Total reaction force at each fixed node.
    """
    boundary_conditions = get_boundary_conditions(nodes)
    reaction_forces = {node.id: np.array([0.0, 0.0, 0.0]) for node in nodes if node.id in boundary_conditions}
    
    for element in elements:
        force = element.internal_force
        direction = np.array([
            (element.node2.x - element.node1.x),
            (element.node2.y - element.node1.y),
            (element.node2.z - element.node1.z)
        ])
        length = np.linalg.norm(direction)
        direction /= length  # Normalize direction vector
        
        force_vector = force * direction
        
        if element.node1.id in reaction_forces:
            reaction_forces[element.node1.id] += force_vector
        if element.node2.id in reaction_forces:
            reaction_forces[element.node2.id] -= force_vector
    
    return reaction_forces


def sum_reaction_forces(reaction_forces):
    """
    Sum total reaction forces across all fixed nodes.
    
    Parameters:
        reaction_forces (dict): Dictionary of reaction forces at fixed nodes.
    
    Returns:
        np.array: Total reaction force vector [Fx, Fy, Fz].
    """
    total_force = np.sum(list(reaction_forces.values()), axis=0)
    return total_force