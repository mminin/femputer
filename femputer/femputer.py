import numpy as np


# Define Node class
class Node:
    def __init__(self, id, x, y):
        self.id = id
        self.x = x
        self.y = y
        self.force = np.array([0.0, 0.0])
        self.displacement = np.array([None, None])
        self.fixed = {"x": False, "y": False}
        
    @staticmethod
    def create_nodes(data): 
        # Initialize nodes from JSON data
        nodes = []
        for id, coords in data["nodes"].items():
            nodes.append(Node(id, coords["x"], coords["y"]))
        # Apply boundary conditions from JSON data
        for node_id, dofs in data["boundary_conditions"].items():
            node = next(node for node in nodes if node.id == node_id)
            for dof in dofs:
                node.fixed[dof] = True
        # Apply forces from JSON data
        for node_id, forces in data["forces"].items():
            node = next(node for node in nodes if node.id == node_id)
            node.force = np.array([forces.get("x", 0.0), forces.get("y", 0.0)])
        return nodes

# Define Element class
class Element:
    def __init__(self, id, node1, node2, material):
        self.id = id
        self.node1 = node1
        self.node2 = node2
        self.material = material
        self.length = np.sqrt((node2.x - node1.x) ** 2 + (node2.y - node1.y) ** 2)
        self.stiffness_matrix = self.calculate_stiffness()
    
    def calculate_stiffness(self):
        A = self.material.area
        E = self.material.young_modulus
        L = self.length
        k = (A * E) / L
        cos_theta = (self.node2.x - self.node1.x) / L
        sin_theta = (self.node2.y - self.node1.y) / L
        k_local = k * np.array(
            [
                [cos_theta ** 2, cos_theta * sin_theta, -cos_theta ** 2, -cos_theta * sin_theta],
                [cos_theta * sin_theta, sin_theta ** 2, -cos_theta * sin_theta, -sin_theta ** 2],
                [-cos_theta ** 2, -cos_theta * sin_theta, cos_theta ** 2, cos_theta * sin_theta],
                [-cos_theta * sin_theta, -sin_theta ** 2, cos_theta * sin_theta, sin_theta ** 2]
            ]
        )
        return k_local

    @staticmethod
    def create_elements(data, nodes, materials):
        # Initialize elements from JSON data
        elements = []
        for id, element_data in data["elements"].items():
            node1 = next(node for node in nodes if node.id == element_data["nodes"][0])
            node2 = next(node for node in nodes if node.id == element_data["nodes"][1])
            material = materials[element_data["material"]]
            elements.append(Element(id, node1, node2, material))
        return elements


# Define Material class
class Material:
    def __init__(self, young_modulus, area):
        self.young_modulus = young_modulus
        self.area = area
        
    @staticmethod
    def create_materials(data):
        materials = {m_id: Material(**props) for m_id, props in data["materials"].items()}
        return materials


