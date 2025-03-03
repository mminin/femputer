import json

def export_to_json(nodes, elements, filename, resolution, forces = None, boundary_conditions=None, intrusion_nodes=None):
    """Export grid data to a JSON file with the updated boundary conditions."""
    data = {
        "nodes": {str(i): {"x": node.x, "y": node.y} for i, node in enumerate(nodes)},
        "elements": {
            str(i): {
                "nodes": [str(element.node1.id), str(element.node2.id)],
                "shape": "default_shape",
                "material": element.material.key  # Use the material attribute of each Element
            }
            for i, element in enumerate(elements)
        },
        "shapes": {
            "default_shape": {
                "area": 0.0001,
                "moment_of_inertia": 1e-06
            }
        },
        "materials": {
            "default_material": {
                "young_modulus": 200000000000.0,
                "tensile_strength": 250000000,
                "compressive_strength": 150000000
            },
            "intrusion_material": {
                "young_modulus": 400000000000.0,
                "tensile_strength": 500000,
                "compressive_strength": 300000000
            }
        },
        "boundary_conditions": boundary_conditions,
        "forces": forces
    }

    with open(filename, 'w') as f:
        json.dump(data, f, indent=4)
    print(f"Grid data exported to {filename}")
