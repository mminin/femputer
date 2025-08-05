import numpy as np
from shapely.geometry import Point, Polygon
from femputer.femputer import Node, Element

def create_grid(x_range, y_range, resolution):
    """Create a 2D grid of nodes."""
    x = np.linspace(x_range[0], x_range[1], resolution[0])
    y = np.linspace(y_range[0], y_range[1], resolution[1])
    xv, yv = np.meshgrid(x, y)
    nodes = [Node(ID=i, x=float(x), y=float(y)) for i, (x, y) in enumerate(zip(xv.ravel(), yv.ravel()))]
    return nodes
    

def select_core_nodes_in_polygon(nodes, intrusion_polygon, radius=10):
    """
    Find nodes inside the intrusion region and select a core area to fix.
    Args:
    - nodes: List of Node objects.
    - intrusion_polygon: Shapely Polygon object representing the intrusion.
    - radius: The radius (in terms of nodes) around the intrusion center to be fixed.
    
    Returns:
    - core_nodes: List of Node objects within the intrusion core.
    """
    # Calculate centroid of the polygon
    centroid = intrusion_polygon.centroid

    # Select nodes within the intrusion area
    core_nodes = []
    for node in nodes:
        if intrusion_polygon.contains(Point(node.x, node.y)):
            # Optional: Apply a distance condition (radius) around the centroid to define the core region
            distance = np.sqrt((node.x - centroid.x) ** 2 + (node.y - centroid.y) ** 2)
            if distance <= radius:  # Core region within radius of centroid
                core_nodes.append(node)
                
    return core_nodes
    
def select_nodes_in_polygon(nodes, polygon):
    selected_nodes = []
    for node in nodes:
        point = Point(node.x, node.y)
        if polygon.contains(point):
            selected_nodes.append(node)
    return selected_nodes
    

def create_elements(nodes, resolution, default_shape, default_material):
    elements = set()
    # resolution = (nx, ny)  # So reshape as (ny, nx) for meshgrid's order
    nx, ny = resolution
    node_indices = np.arange(len(nodes)).reshape((ny, nx))
    for i in range(ny - 1):     # rows (y)
        for j in range(nx - 1): # columns (x)
            n1 = node_indices[i, j]
            n2 = node_indices[i, j + 1]
            n3 = node_indices[i + 1, j + 1]
            n4 = node_indices[i + 1, j]
            
            elements.add((n1, n2))  # Bottom edge (horizontal)
            elements.add((n2, n3))  # Right edge   (vertical)
            elements.add((n3, n4))  # Top edge     (horizontal)
            elements.add((n4, n1))  # Left edge    (vertical)
            elements.add((n1, n3))  # Diagonal
            elements.add((n2, n4))  # Diagonal

    element_objects = [
        Element(ID=i, node1=nodes[n1], node2=nodes[n2], shape=default_shape, material=default_material)
        for i, (n1, n2) in enumerate(elements)
    ]
    return element_objects
