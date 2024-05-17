import bpy



def initialize_material(name):
    material = bpy.data.materials.new(name=name)
    material.use_nodes = True
    nodes = material.node_tree.nodes
    links = material.node_tree.links
    nodes.clear()
    return material, nodes, links


def create_node_group(name='MyCustomNode', start=(0, 0), end=(300, 0)):
    # Create a new node group
    node_group = bpy.data.node_groups.new(type="ShaderNodeTree", name=name)

    # Nodes are added to the node group's own nodes collection
    group_inputs = node_group.nodes.new("NodeGroupInput")
    group_inputs.location = start

    group_outputs = node_group.nodes.new("NodeGroupOutput")
    group_outputs.location = end
    
    return node_group, group_inputs, group_outputs


def add_node_group(nodes, node_group, pos=(0, 0)):
    group_node = nodes.new("ShaderNodeGroup")
    group_node.node_tree = node_group
    group_node.location = pos
    return group_node


def add_node_group_name(nodes, name, pos=(0, 0)):
    group_node = nodes.new("ShaderNodeGroup")
    group_node.node_tree = bpy.data.node_groups[name]
    group_node.location = pos
    return group_node
