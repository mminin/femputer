import numpy as np
import pandas as pd
from pathlib import Path


# Define Node class
class Node:
    def __init__(self, ID, x, y, z=0.0):
        self.id = ID
        self.x = x
        self.y = y
        self.z = z
        self.force = np.array([0.0, 0.0, 0.0])
        self.displacement = np.array([None, None, None])
        self.fixed = {"x": False, "y": False, "z": False}
    
    @staticmethod
    def create_nodes(data, is_3d=True):
        """Initialize nodes from JSON data in 2D or 3D."""
        nodes = []
        for ID, coords in data["nodes"].items():
            nodes.append(Node(ID, coords["x"], coords["y"], coords.get("z", 0.0) if is_3d else 0.0))
        
        # Apply boundary conditions
        for node_id, dofs in data["boundary_conditions"].items():
            node = next(node for node in nodes if node.id == node_id)
            for dof in dofs:
                node.fixed[dof] = True
        
        # Apply forces
        for node_id, forces in data["forces"].items():
            node = next(node for node in nodes if node.id == node_id)
            node.force = np.array([
                forces.get("x", 0.0),
                forces.get("y", 0.0),
                forces.get("z", 0.0) if is_3d else 0.0
            ])
        return nodes

# Define Element class
class Element:
    def __init__(self, ID, node1, node2, shape, material):
        self.id = ID
        self.node1 = node1
        self.node2 = node2
        self.shape = shape
        self.material = material
        self.length = np.sqrt((node2.x - node1.x) ** 2 + (node2.y - node1.y) ** 2 + (node2.z - node1.z) ** 2)
        self.stiffness_matrix = self.calculate_stiffness()
        
    def update_properties(self):
        """Update element length and stiffness matrix."""
        #self.length = ((self.node2.x - self.node1.x)**2 + (self.node2.y - self.node1.y)**2)**0.5
        self.length = np.sqrt((self.node2.x - self.node1.x) ** 2 + (self.node2.y - self.node1.y) ** 2 + (self.node2.z - self.node1.z) ** 2)
        self.stiffness_matrix = self.calculate_stiffness()

    def calculate_stiffness(self):
        """Calculate the stiffness matrix for the element."""
        A = self.shape.area  # Access area from the Shape object
        E = self.material.young_modulus  # Access Young's modulus from the Material object
        L = self.length
        #print(self.id)
        #print(self.node1.x, self.node1.y, self.node1.z)
        #print(self.node2.x, self.node2.y, self.node2.z)
        #print(L)
        k = (A * E) / L
        cos_theta = (self.node2.x - self.node1.x) / L
        sin_theta = (self.node2.y - self.node1.y) / L
        sin_phi = (self.node2.z - self.node1.z) / L
        k_local = k * np.array([
            [cos_theta ** 2, cos_theta * sin_theta, cos_theta * sin_phi, -cos_theta ** 2, -cos_theta * sin_theta, -cos_theta * sin_phi],
            [cos_theta * sin_theta, sin_theta ** 2, sin_theta * sin_phi, -cos_theta * sin_theta, -sin_theta ** 2, -sin_theta * sin_phi],
            [cos_theta * sin_phi, sin_theta * sin_phi, sin_phi ** 2, -cos_theta * sin_phi, -sin_theta * sin_phi, -sin_phi ** 2],
            [-cos_theta ** 2, -cos_theta * sin_theta, -cos_theta * sin_phi, cos_theta ** 2, cos_theta * sin_theta, cos_theta * sin_phi],
            [-cos_theta * sin_theta, -sin_theta ** 2, -sin_theta * sin_phi, cos_theta * sin_theta, sin_theta ** 2, sin_theta * sin_phi],
            [-cos_theta * sin_phi, -sin_theta * sin_phi, -sin_phi ** 2, cos_theta * sin_phi, sin_theta * sin_phi, sin_phi ** 2]
        ])
        #print("stiffness mtx, \n")
        #print(k_local.round(2))
        return k_local
    
    def stiffness_matrix_3d(self):
        """Alias for calculate_stiffness to match function call in femsolver.py."""
        return self.calculate_stiffness()
    
    def stiffness_matrix_2d(self):
        """Compute a 2D version of the stiffness matrix by removing z-components."""
        full_matrix = self.calculate_stiffness()
        return full_matrix[np.ix_([0,1,3,4],[0,1,3,4])]


    @staticmethod
    def set_material(material_id, materials):
        """
        Retrieve the material for an element.

        Parameters:
            material_id (str): ID of the material.
            materials (dict): Dictionary of predefined materials.

        Returns:
            Material: A Material object with the appropriate properties.

        Raises:
            ValueError: If the material ID is undefined.
        """

        # Case 1: Material is predefined in the JSON data
        if material_id in materials:
            return materials[material_id]

        # Case 2: Material ID is undefined
        raise ValueError(f"Material '{material_id}' is undefined. Available materials: {list(materials.keys())}")


    @staticmethod
    def set_shape(shape_id, shapes, lookup_table):
        """
        Retrieve or define the shape for an element.

        Parameters:
            shape_id (str): ID of the shape.
            shapes (dict): Dictionary of predefined shapes.
            lookup_table (DataFrame): DataFrame containing shape properties.

        Returns:
            Shape: A Shape object with the appropriate properties.

        Raises:
            ValueError: If the shape ID is undefined.
        """

        # Case 1: Shape is predefined in the JSON data
        if shape_id in shapes:
            return shapes[shape_id]

        # Case 2: Attempt to lookup shape in the external table
        shape_row = lookup_table.loc[lookup_table["Shape"] == shape_id]
        if not shape_row.empty:
            return Shape(
                key=shape_id,
                area=shape_row["Area (m2)"].values[0],
                moment_of_inertia=shape_row.get("Ix (m4)", None).values[0]
            )

        # Case 3: Shape ID is undefined
        raise ValueError(f"Shape '{shape_id}' is undefined. "
                        f"Available shapes: {list(shapes.keys()) + lookup_table['Shape'].tolist()}")

    @staticmethod
    def create_elements(data, nodes, shapes, materials, lt_fname):
        """
        Create elements by assigning nodes, shapes, and materials.

        Parameters:
            data (dict): JSON data containing element definitions.
            nodes (list): List of Node objects.
            shapes (dict): Dictionary of Shape objects.
            materials (dict): Dictionary of Material objects.

        Returns:
            list: A list of Element objects.
        """
        # Load the lookup table for shapes
        script_dir = Path(__file__).resolve().parent
        #lookup_table_path = script_dir / "LookupTables" / "beam_pricing_all.csv"
        lookup_table_path = lt_fname
        lookup_table = pd.read_csv(lookup_table_path)

        elements = []
        for ID, element_data in data["elements"].items():
            # Assign nodes
            node1 = next(node for node in nodes if node.id == element_data["nodes"][0])
            node2 = next(node for node in nodes if node.id == element_data["nodes"][1])

            # Assign shape using the new `set_shape` method
            shape = Element.set_shape(element_data["shape"], shapes, lookup_table)

            # Assign material using the updated `set_material` method
            material = Element.set_material(element_data["material"], materials)

            # Create and append the element
            elements.append(Element(ID, node1, node2, shape, material))

        return elements


class Shape:
    def __init__(self, key, area, moment_of_inertia=None): # Todo: Provide Ix and Iy instead of MoI
        self.key = key
        self.area = area
        self.moment_of_inertia = moment_of_inertia

    @staticmethod
    def create_shapes(data):
        desired_keys = {"area", "moment_of_inertia"}  # List desired keys
        shapes = {
            s_id: Shape(
                key=s_id,
                **{k: v for k, v in props.items() if k in desired_keys}  # Filter for desired keys
            )
            for s_id, props in data["shapes"].items()
        }
        return shapes

class Material:
    def __init__(self, key, young_modulus, tensile_strength, compressive_strength):
        self.key = key
        self.young_modulus = young_modulus
        self.tensile_strength = tensile_strength
        self.compressive_strength = compressive_strength

    @staticmethod
    def create_materials(data):
        desired_keys = {"young_modulus", "tensile_strength", "compressive_strength"}  # List desired keys
        materials = {
            m_id: Material(
                key=m_id,
                **{k: v for k, v in props.items() if k in desired_keys}  # Filter for desired keys
            )
            for m_id, props in data["materials"].items()
        }
        return materials

def create_and_assign_material(element, selected_beam, data):
    """
    Create a new shape for the selected beam and assign it to the element.

    Parameters:
        element: The element to which the shape will be assigned.
        selected_beam: The selected beam data (row from beam_table).
        data: The JSON data for adding new shapes.
    """
    # Create a new shape for the selected beam
    new_shape_key = selected_beam['Shape']
    new_shape = Shape(
        key=new_shape_key,
        area=selected_beam['Area (m2)'],  # Update to the selected beam's area
        moment_of_inertia=selected_beam.get('Ix (m4)', 0.0)  # Default moment of inertia if not provided
    )
    
    # Add the new shape to the JSON data for export
    if "shapes" not in data:
        data["shapes"] = {}
    data["shapes"][new_shape_key] = {
        "area": new_shape.area,
        "moment_of_inertia": new_shape.moment_of_inertia
    }

    # Update the element's shape properties
    element.shape = new_shape

    return new_shape_key

def select_beam(elements, beam_table, data, factor_of_safety=4):
    K=1
    beam_assignments = {}
    for element in elements:
        F = abs(element.internal_force)  # Internal force (N)
        L = element.length  # Element length (m)
        E = element.material.young_modulus  # Young's modulus (Pa)
        material_yield_strength = 250000000  # Yield strength for A36 steel

        # Yielding check
        A_req_yield = F / (material_yield_strength / factor_of_safety)
        print(F, material_yield_strength, factor_of_safety, round(A_req_yield, 7))
        valid_beams = beam_table[beam_table['Area (m2)'] >= A_req_yield].copy()
        print("Yielding only solution is ", valid_beams.sort_values("Cost per Meter (usd)").iloc[0])

        if valid_beams.empty:
            raise RuntimeError(f"No valid beams for element {element.id} under yielding criteria.")

        if element.internal_force < 0:  # Compressive force -> check for buckling
            # Buckling check, use minor radius of gyration, Iy - as this is the one that will buckle first
            valid_beams['Radius of Gyration (m)'] = (valid_beams['Iy (m4)'] / valid_beams['Area (m2)'])**0.5
            slenderness_ratio = (K * L) / valid_beams['Radius of Gyration (m)']

            # Critical slenderness ratio
            critical_slenderness = (np.pi**2 * E / material_yield_strength)**0.5

            # Determine critical stress based on slenderness ratio
            valid_beams['Critical Stress (Pa)'] = np.where(
                slenderness_ratio <= critical_slenderness,
                # Stocky columns (Johnson's formula)
                material_yield_strength * (1 - (material_yield_strength / (4 * (np.pi**2 * E)))),
                # Slender columns (Euler's formula)
                (np.pi**2 * E) / (slenderness_ratio**2)
            )

            # Filter beams that meet buckling criteria
            valid_beams = valid_beams[valid_beams['Critical Stress (Pa)'] >= (F / valid_beams['Area (m2)'])]
            print("Buckling possible, new solution is ", valid_beams.sort_values("Cost per Meter (usd)").iloc[0])

        if valid_beams.empty:
            raise RuntimeError(f"No valid beams for element {element.id} under buckling criteria.")

        # Select cheapest beam
        selected_beam = valid_beams.sort_values("Cost per Meter (usd)").iloc[0]
        # Delegate the creation and assignment of the material
        beam_assignments[element.id] = create_and_assign_material(element, selected_beam, data)

    return beam_assignments
