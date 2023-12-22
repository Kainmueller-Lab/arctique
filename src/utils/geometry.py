from scipy.spatial import Voronoi # NOTE: You might install it directly into Blender's python using pip. - ck


def print_voronoi_stats(vor):    
    # Vertex locations
    print(f"Vertices: {vor.vertices}")
    # - The Voronoi vertices represent the points where three or more Voronoi cells meet.
    # - Each vertex is a 3D point in the input space.

    # Ridge vertices
    print(f"Ridge vertices: {vor.ridge_vertices}")
    # - The ridge vertices represent the vertices of the Voronoi ridges (faces between Voronoi cells).
    # - Each ridge is represented by a set of vertex indices.
    # - If one of the indices is -1, the ridge extends to infinity.

    # Ridge dictionary
    print(f"Ridge dict: {vor.ridge_dict}")
    # - The ridge dictionary maps pairs of point indices to the corresponding ridge index.
    # - It provides a way to look up the ridge index based on the indices of the points it connects.

    # Points
    print(f"Regions per point: {vor.point_region}")
    # - The point_region array maps each input point to the index of the Voronoi region (cell) it belongs to.
    # - The index corresponds to the index of the Voronoi region in the 'regions' array.

    # Ridge points
    print(f"Ridge points: {vor.ridge_points}")
    # - The ridge_points array provides the indices of the points that form the corners of each Voronoi ridge.
    # - Each row corresponds to a ridge, and the two columns contain the indices of the points.

    # Regions
    print(f"Regions: {vor.regions}")
    # - The regions array contains the vertices of each Voronoi region (polygonal cell).
    # - It includes the indices of the Voronoi vertices that form the boundary of the region.
    # - Some entries may be empty lists, representing unbounded regions.

def finite_region_points(vor):
    fr_points = []
    for point_idx, region_idx in enumerate(vor.point_region):
        region = vor.regions[region_idx]
        if -1 not in region and len(region) > 0:
            fr_points.append(point_idx)
    return fr_points

def compute_faces_by_seeds(vor, seeds):
    ridges = {}
    # Print ridge pairs of finite regions
    for point_idx in seeds:
        # Get boundary ridges (= faces, as list of face vertices)
        ridge_list = []
        for ridge_pair in vor.ridge_dict.keys():
            if point_idx in set(ridge_pair):
                ridge = vor.ridge_dict[ridge_pair]
                #print(f"Boundary ridge pair: {ridge_pair} - {ridge}")
                ridge_list.append(ridge)
        ridges[point_idx] = ridge_list
    print(f"Created the Voronoi region faces")
    return ridges

def voronoi_test():
    # Generate 27 integer points in the 0, -1, 1 lattice
    points = [[i, j, k] for i in range(-1, 2) for j in range(-1, 2) for k in range(-1, 2)]
    vor = Voronoi(points)
    fr_points = finite_region_points(vor)
    assert len(fr_points) == 1
    assert fr_points[0] == 13