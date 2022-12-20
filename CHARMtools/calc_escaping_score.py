from email import parser
import open3d as o3d
import numpy as np
import pandas as pd
import argparse
import copy

parser = argparse.ArgumentParser()
parser.add_argument('-i','--input', type=str, help='input point cloud 3dg file')
parser.add_argument('-a','--alpha', type=str, help='parameter alpha for alpha shape',default=2)
parser.add_argument('-o','--output', type=str, help='output distance to lamin(alpha shape")')
args = parser.parse_args()

# Load point cloud
test_structure = pd.read_csv(args.input ,sep="\t",header=None)

alpha = args.alpha

pcd = o3d.geometry.PointCloud()
pcd.points = o3d.utility.Vector3dVector(test_structure[[2,3,4]].values)

print(f"alpha={alpha:.3f}")
mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_alpha_shape(pcd, alpha)
mesh.compute_vertex_normals()

print("Cluster connected triangles")
with o3d.utility.VerbosityContextManager(
        o3d.utility.VerbosityLevel.Debug) as cm:
    triangle_clusters, cluster_n_triangles, cluster_area = (
        mesh.cluster_connected_triangles())
triangle_clusters = np.asarray(triangle_clusters)
cluster_n_triangles = np.asarray(cluster_n_triangles)
cluster_area = np.asarray(cluster_area)

print("Meesh with small clusters (<500) will be removed")
mesh_0 = copy.deepcopy(mesh)
triangles_to_remove = cluster_n_triangles[triangle_clusters] < 500
mesh_0.remove_triangles_by_mask(triangles_to_remove)

mesh_0.compute_vertex_normals()

meshn = o3d.t.geometry.TriangleMesh.from_legacy(mesh_0)

# Create a scene and add the triangle mesh
scene = o3d.t.geometry.RaycastingScene()
_ = scene.add_triangles(meshn)  # we do not need the geometry ID for mesh

query_point = o3d.core.Tensor([test_structure[[2,3,4]].values.tolist()], dtype=o3d.core.Dtype.Float32)

# Compute distance of the query point from the surface
unsigned_distance = scene.compute_distance(query_point)
signed_distance = scene.compute_signed_distance(query_point)
occupancy = scene.compute_occupancy(query_point)

test_structure["distance"] = unsigned_distance[0].numpy()
test_structure.columns = ["chrom","pos","x","y","z","distance"]
test_structure[["chrom","pos","distance"]].to_csv(args.output,sep="\t",index=False)