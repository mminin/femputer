from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
import numpy as np
import json
import pandas as pd
from .femputer import Node, Element, Material, Shape

def perform_fem_analysis(nodes, elements):
    """Run FEM analysis: stiffness matrix assembly, boundary conditions, and stress computation."""
    # Assemble global stiffness matrix
    K = assemble_stiffness_matrix(nodes, elements)
    apply_boundary_conditions_and_solve(K, nodes)
    compute_stresses(elements)
    return nodes, elements

class Solver:
    def __init__(self, json_file, beam_csv, pricing_csv):
        self.json_file = json_file
        self.beam_csv = beam_csv
        self.pricing_csv = pricing_csv
        self.data = None
        self.beam_table = None
        self.pricing_table = None

    def load_input_data(self, json_file=None, beam_csv=None, pricing_csv=None):
        """Load input data for nodes, elements, shapes, and materials."""
        json_file = json_file or self.json_file
        beam_csv = beam_csv or self.beam_csv
        pricing_csv = pricing_csv or self.pricing_csv

        # Load JSON data
        with open(json_file, 'r') as f:
            self.data = json.load(f)
        
        # Load and standardize beam table
        self.beam_table = pd.read_csv(beam_csv)

        # Load pricing table
        self.pricing_table = pd.read_csv(pricing_csv)
        return self.data, self.beam_table, self.pricing_table

    def preprocess_data(self, data=None):
        """Initialize nodes, shapes, materials, and elements."""
        data = data or self.data  # Default to the loaded data
        if data is None:
            raise ValueError("No data provided for preprocessing.")

        # Initialize nodes, shapes, and materials
        nodes = Node.create_nodes(data)
        shapes = Shape.create_shapes(data)  # Load shapes from JSON
        materials = Material.create_materials(data)  # Load materials from JSON

        # Initialize elements with nodes, shapes, and materials
        elements = Element.create_elements(data, nodes, shapes, materials)
        return nodes, elements

def assemble_stiffness_matrix(nodes, elements):
    num_dofs = len(nodes) * 2
    K = lil_matrix((num_dofs, num_dofs))
    
    for element in elements:
        n1 = nodes.index(element.node1)
        n2 = nodes.index(element.node2)
        k_local = element.stiffness_matrix
        
        dof_indices = [n1 * 2, n1 * 2 + 1, n2 * 2, n2 * 2 + 1]
        for i in range(4):
            for j in range(4):
                K[dof_indices[i], dof_indices[j]] += k_local[i, j]
    
    return K

def apply_boundary_conditions_and_solve(K, nodes):
    num_dofs = len(nodes) * 2
    F = np.zeros(num_dofs)
    fixed_dofs = []
    
    # Apply loads and boundary conditions
    for node in nodes:
        idx = nodes.index(node) * 2
        if node.fixed["x"]:
            fixed_dofs.append(idx)
        else:
            F[idx] = node.force[0]
        if node.fixed["y"]:
            fixed_dofs.append(idx + 1)
        else:
            F[idx + 1] = node.force[1]
    
    # Modify stiffness matrix for fixed DOFs
    for dof in fixed_dofs:
        K[dof, :] = 0
        K[:, dof] = 0
        K[dof, dof] = 1
        F[dof] = 0
    
    # Solve for displacements
    displacements = spsolve(K.tocsr(), F)
    
    # Assign displacements to nodes
    for node in nodes:
        idx = nodes.index(node) * 2
        node.displacement = np.array([displacements[idx], displacements[idx + 1]])
        
def compute_stresses(elements):
    """
    Compute stresses and internal forces for each element.

    Parameters:
        elements (list): List of Element objects.
    """
    for element in elements:
        # Calculate strain
        length = element.length
        displacement_diff = element.node2.displacement - element.node1.displacement
        cos_theta = (element.node2.x - element.node1.x) / length
        sin_theta = (element.node2.y - element.node1.y) / length
        axial_disp = displacement_diff[0] * cos_theta + displacement_diff[1] * sin_theta
        strain = axial_disp / length

        # Calculate stress
        stress = element.material.young_modulus * strain
        element.stress = stress

        # Calculate internal force
        element.internal_force = element.stress * element.shape.area


def convert_to_metric(json_file, output_file):
    """Convert imperial units in the input JSON file to metric units."""
    with open(json_file, 'r') as f:
        data = json.load(f)

    if data.get("units", "metric") == "imperial":
        # Conversion factors
        inch_to_m = 0.0254
        feet_to_m = 0.3048
        psi_to_pa = 6894.76
        lbf_to_n = 4.44822
        in2_to_m2 = 0.00064516

        # Convert node coordinates
        for node in data["nodes"].values():
#            node["x"] *= inch_to_m
#            node["y"] *= inch_to_m
            node["x"] *= feet_to_m
            node["y"] *= feet_to_m

        # Convert forces
        for force in data["forces"].values():
            if "x" in force:
                force["x"] *= lbf_to_n
            if "y" in force:
                force["y"] *= lbf_to_n

        # Convert material properties
        for material in data["materials"].values():
            material["young_modulus"] *= psi_to_pa
            material["area"] *= in2_to_m2

        data["units"] = "metric"
        with open(output_file, 'w') as out_f:
            json.dump(data, out_f, indent=4)
        print(f"Converted data saved to {output_file}.")
    else:
        print("Units are already in metric. No conversion needed.")



