from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
import numpy as np

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
    for element in elements:
        length = element.length
        displacement_diff = element.node2.displacement - element.node1.displacement
        cos_theta = (element.node2.x - element.node1.x) / length
        sin_theta = (element.node2.y - element.node1.y) / length
        axial_disp = displacement_diff[0] * cos_theta + displacement_diff[1] * sin_theta
        strain = axial_disp / length
        stress = element.material.young_modulus * strain
        element.stress = stress
        #print(f"Element {element.id}: Stress = {element.stress} Pa")
