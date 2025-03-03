from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
import numpy as np
import json
import pandas as pd
from .femputer import Node, Element, Material, Shape

def perform_fem_analysis(nodes, elements, is_3d=True):
    """Run FEM analysis in 2D or 3D: stiffness matrix assembly, boundary conditions, and stress computation."""
    K = assemble_stiffness_matrix(nodes, elements, is_3d)
    apply_boundary_conditions_and_solve(K, nodes, is_3d)
    compute_stresses(elements, is_3d)
    return nodes, elements

class Solver:
    def __init__(self, json_file, beam_csv, pricing_csv):
        self.json_file = json_file
        self.beam_csv = beam_csv
        self.pricing_csv = pricing_csv
        self.data = None
        self.beam_table = None
        self.pricing_table = None

    def load_input_data(self):
        """Load input data from JSON and CSV files."""
        with open(self.json_file, 'r') as f:
            self.data = json.load(f)
        self.beam_table = pd.read_csv(self.beam_csv)
        self.pricing_table = pd.read_csv(self.pricing_csv)
        return self.data, self.beam_table, self.pricing_table

    def preprocess_data(self):
        """Initialize nodes, elements, shapes, and materials from data."""
        if self.data is None:
            raise ValueError("No data loaded.")
        lt_fname = self.pricing_csv
        nodes = Node.create_nodes(self.data)
        shapes = Shape.create_shapes(self.data)
        materials = Material.create_materials(self.data)
        elements = Element.create_elements(self.data, nodes, shapes, materials, lt_fname)
        return nodes, elements

def assemble_stiffness_matrix(nodes, elements, is_3d=True):
    num_dofs_per_node = 3 if is_3d else 2
    num_dofs = len(nodes) * num_dofs_per_node
    K = lil_matrix((num_dofs, num_dofs))
    
    for element in elements:
        n1 = nodes.index(element.node1)
        n2 = nodes.index(element.node2)
        k_local = element.stiffness_matrix_3d() if is_3d else element.stiffness_matrix_2d()
        
        dof_indices = []
        for n in [n1, n2]:
            dof_indices.extend([n * num_dofs_per_node + i for i in range(num_dofs_per_node)])
        
        for i, global_i in enumerate(dof_indices):
            for j, global_j in enumerate(dof_indices):
                K[global_i, global_j] += k_local[i, j]
    
    return K

def apply_boundary_conditions_and_solve(K, nodes, is_3d=True):
    num_dofs_per_node = 3 if is_3d else 2
    num_dofs = len(nodes) * num_dofs_per_node
    F = np.zeros(num_dofs)
    fixed_dofs = []
    
    for node in nodes:
        base_idx = nodes.index(node) * num_dofs_per_node
        for i, axis in enumerate(["x", "y", "z"][:num_dofs_per_node]):
            if node.fixed[axis]:
                fixed_dofs.append(base_idx + i)
            else:
                F[base_idx + i] = node.force[i]
    
    for dof in fixed_dofs:
        K[dof, :] = 0
        K[:, dof] = 0
        K[dof, dof] = 1
        F[dof] = 0
    
    displacements = spsolve(K.tocsr(), F)
    
    for node in nodes:
        base_idx = nodes.index(node) * num_dofs_per_node
        node.displacement = np.array([displacements[base_idx + i] for i in range(num_dofs_per_node)])

def compute_stresses(elements, is_3d=True):
    num_dims = 3 if is_3d else 2
    for element in elements:
        length = element.length
        displacement_diff = element.node2.displacement - element.node1.displacement
        direction_vector = np.array([
            (element.node2.x - element.node1.x),
            (element.node2.y - element.node1.y),
            *( [(element.node2.z - element.node1.z)] if is_3d else [] )
        ]) / length
        
        axial_disp = np.dot(displacement_diff, direction_vector)
        strain = axial_disp / length
        stress = element.material.young_modulus * strain
        element.stress = stress
        element.internal_force = stress * element.shape.area


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
        in4_to_m4 = (0.0254 ** 4)  # Convert in^4 to m^4 for moments of inertia


        # Convert node coordinates
        for node in data["nodes"].values():
#            node["x"] *= inch_to_m
#            node["y"] *= inch_to_m
            node["x"] *= feet_to_m
            node["y"] *= feet_to_m
            if "z" in node:
                node["z"] *= inch_to_m

        # Convert forces
        for force in data["forces"].values():
            if "x" in force:
                force["x"] *= lbf_to_n
            if "y" in force:
                force["y"] *= lbf_to_n
            if "z" in force:
                force["z"] *= lbf_to_n

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



