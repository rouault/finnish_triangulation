from osgeo import ogr, osr
import json
import random
import os

# Customization of behaviour
use_refinement = True


dir_data = os.path.dirname(__file__)

kkjproj = osr.SpatialReference()
if hasattr(kkjproj, 'SetAxisMappingStrategy'):
    kkjproj.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
kkjproj.ImportFromEPSG(2393)
kkjgeog = osr.SpatialReference()
if hasattr(kkjgeog, 'SetAxisMappingStrategy'):
    kkjgeog.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
kkjgeog.ImportFromEPSG(4123)
kkjproj_to_kkjgeog = osr.CoordinateTransformation(kkjproj, kkjgeog)
kkjgeog_to_kkjproj = osr.CoordinateTransformation(kkjgeog, kkjproj)

etrs35fin = osr.SpatialReference()
if hasattr(etrs35fin, 'SetAxisMappingStrategy'):
    etrs35fin.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
etrs35fin.ImportFromEPSG(3067)
etrsgeog = osr.SpatialReference()
if hasattr(etrsgeog, 'SetAxisMappingStrategy'):
    etrsgeog.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
etrsgeog.ImportFromEPSG(4258)
etrs35fin_to_etrsgeog = osr.CoordinateTransformation(etrs35fin, etrsgeog)
etrsgeog_to_etrs35fin = osr.CoordinateTransformation(etrsgeog, etrs35fin)

def get_barycentric_coord(x1, y1, x2, y2, x3, y3, x, y):
    det_T = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)
    lambda_1 = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / det_T
    lambda_2 = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / det_T
    lambda_3 = 1 - lambda_1 - lambda_2
    return lambda_1, lambda_2, lambda_3

def find_triangle(triangles, vertices, x, y):
    for (id1, id2, id3) in triangles:
        x1 = vertices[id1][0]
        y1 = vertices[id1][1]
        x2 = vertices[id2][0]
        y2 = vertices[id2][1]
        x3 = vertices[id3][0]
        y3 = vertices[id3][1]
        lambda_1, lambda_2, lambda_3 = get_barycentric_coord(x1, y1, x2, y2, x3, y3, x, y)
        if lambda_1 >= -1e-10 and lambda_1 <= 1+1e-10 and \
           lambda_2 >= -1e-10 and lambda_2 <= 1+1e-10 and \
           lambda_3 >= -1e-10 and lambda_3 <= 1+1e-10:
               return id1, id2, id3, lambda_1, lambda_2, lambda_3
    return None

def read_vertices():
    vertices = {}
    vertices_geog = {}
    with open(os.path.join(dir_data, 'kkjEUREFFINtriangulationVertices.txt'), 'rt') as f:
        for line in f.readlines():
            tri_id, ykj_north, ykj_east, tm35fin_north, tm35fin_east = line.rstrip().split('\t')
            ykj_north = float(ykj_north)
            ykj_east = float(ykj_east)
            tm35fin_north = float(tm35fin_north)
            tm35fin_east = float(tm35fin_east)
            vertices[tri_id] = [ykj_east, ykj_north, tm35fin_east, tm35fin_north]

            kkj_lon, kkj_lat, _ = kkjproj_to_kkjgeog.TransformPoint(ykj_east, ykj_north)
            etrs_lon, etrs_lat, _ = etrs35fin_to_etrsgeog.TransformPoint(tm35fin_east, tm35fin_north)
            vertices_geog[tri_id] = [kkj_lon, kkj_lat, etrs_lon, etrs_lat]

    return vertices, vertices_geog

def read_triangles():
    triangles = []
    geomcoll = ogr.Geometry(ogr.wkbMultiPolygon)
    with open(os.path.join(dir_data, 'kkjEUREFFINtriangulationNetwork.txt'), 'rt') as f:
        for line in f.readlines():
            id1, id2, id3 = line.rstrip().split(' ')
            assert id1 in vertices
            assert id2 in vertices
            assert id3 in vertices
            triangles.append((id1, id2, id3))
            x1 = vertices[id1][0]
            y1 = vertices[id1][1]
            x2 = vertices[id2][0]
            y2 = vertices[id2][1]
            x3 = vertices[id3][0]
            y3 = vertices[id3][1]
            tri_geom = ogr.CreateGeometryFromWkt('POLYGON((%.18g %.18g,%.18g %.18g,%.18g %.18g,%.18g %.18g))' % (x1,y1,x2,y2,x3,y3,x1,y1))
            geomcoll.AddGeometry(tri_geom)

    envelope = geomcoll.UnionCascaded()
    return triangles, envelope

def read_grid(filename):
    # Header 2393 2393 3067 6600000.0 7800000.0 3050000.0 3750000.0 1000.0 1000.0 1
    with open(filename, 'rt') as f:
        _, _, _, min_northing, max_northing, min_easting, max_easting, step_1, step_2, _ = f.readline().split(' ')
        min_northing = float(min_northing)
        max_northing = float(max_northing)
        min_easting = float(min_easting)
        max_easting = float(max_easting)
        step_1 = float(step_1)
        step_2 = float(step_2)
        assert step_1 == step_2
        lines = f.readlines()
        assert(len(lines)) == 1 + int((max_northing - min_northing) / step_1 + 0.5)
        grid = []
        for j in range(len(lines)):
            vals = lines[len(lines)-1-j].rstrip().split(' ') # read from bottom to top
            assert len(vals) == 1 + int((max_easting - min_easting) / step_1 + 0.5)
            grid.append( [float(v) for v in vals] )
        return min_easting, max_easting, step_1, min_northing, max_northing, step_1, grid


def etrs_exact_from_xy(triangles, vertices, x, y):
    id1, id2, id3, lambda_1, lambda_2, lambda_3 = find_triangle(triangles, vertices, x, y)
    x_tm35fin = vertices[id1][2] * lambda_1 + vertices[id2][2] * lambda_2 + vertices[id3][2] * lambda_3
    y_tm35fin = vertices[id1][3] * lambda_1 + vertices[id2][3] * lambda_2 + vertices[id3][3] * lambda_3
    etrs_lon_exact, etrs_lat_exact, _ = etrs35fin_to_etrsgeog.TransformPoint(x_tm35fin, y_tm35fin)
    return x_tm35fin, y_tm35fin, etrs_lon_exact, etrs_lat_exact

def etrs_exact_from_lonlat(triangles, vertices, lon_kkj, lat_kkj):
    x, y, _ = kkjgeog_to_kkjproj.TransformPoint(lon_kkj, lat_kkj)
    return etrs_exact_from_xy(triangles, vertices, x, y)

def etrs_approx_from_lonlat(vertices_geog, lon1, lat1, lon2, lat2, lon3, lat3, id1_geog, id2_geog, id3_geog, lon_kkj, lat_kkj):
    lambda_1, lambda_2, lambda_3 = get_barycentric_coord(lon1, lat1, lon2, lat2, lon3, lat3, lon_kkj, lat_kkj)
    etrs_lon = vertices_geog[id1_geog][2] * lambda_1 + vertices_geog[id2_geog][2] * lambda_2 + vertices_geog[id3_geog][2] * lambda_3
    etrs_lat = vertices_geog[id1_geog][3] * lambda_1 + vertices_geog[id2_geog][3] * lambda_2 + vertices_geog[id3_geog][3] * lambda_3
    x_tm35fin_from_geog_interp, y_tm35fin_from_geog_interp, _ = etrsgeog_to_etrs35fin.TransformPoint(etrs_lon, etrs_lat)
    return x_tm35fin_from_geog_interp, y_tm35fin_from_geog_interp, etrs_lon, etrs_lat

def refine_triangle(vertices, id1_geog, id2_geog, id3_geog, refinement_threshold_in_shape, shape, refinement_threshold_outside, vertices_geog, new_triangles, map_long_to_short_id, depth=0):
    lon1 = vertices_geog[id1_geog][0]
    lat1 = vertices_geog[id1_geog][1]
    lon2 = vertices_geog[id2_geog][0]
    lat2 = vertices_geog[id2_geog][1]
    lon3 = vertices_geog[id3_geog][0]
    lat3 = vertices_geog[id3_geog][1]

    #print(lon1, lat1, lon2, lat2, lon3, lat3, id1_geog, id2_geog, id3_geog)
    if False:
        x_tm35fin_from_geog_interp, y_tm35fin_from_geog_interp, _, _ = etrs_approx_from_lonlat(vertices_geog, lon1, lat1, lon2, lat2, lon3, lat3, id1_geog, id2_geog, id3_geog, lon1, lat1)
        x_tm35fin, y_tm35fin, _, _ = etrs_exact_from_lonlat(triangles, vertices, lon1, lat1)
        assert abs(x_tm35fin_from_geog_interp - x_tm35fin) < 1e-3, (id1_geog, x_tm35fin_from_geog_interp, x_tm35fin)
        assert abs(y_tm35fin_from_geog_interp - y_tm35fin) < 1e-3, (id1_geog, y_tm35fin_from_geog_interp, y_tm35fin)
        x_tm35fin_from_geog_interp, y_tm35fin_from_geog_interp, _, _ = etrs_approx_from_lonlat(vertices_geog, lon1, lat1, lon2, lat2, lon3, lat3, id1_geog, id2_geog, id3_geog, lon2, lat2)
        x_tm35fin, y_tm35fin, _, _ = etrs_exact_from_lonlat(triangles, vertices, lon2, lat2)
        assert abs(x_tm35fin_from_geog_interp - x_tm35fin) < 1e-3, (id2_geog, x_tm35fin_from_geog_interp, x_tm35fin)
        assert abs(y_tm35fin_from_geog_interp - y_tm35fin) < 1e-3, (id2_geog, y_tm35fin_from_geog_interp, y_tm35fin)
        x_tm35fin_from_geog_interp, y_tm35fin_from_geog_interp, _, _ = etrs_approx_from_lonlat(vertices_geog, lon1, lat1, lon2, lat2, lon3, lat3, id1_geog, id2_geog, id3_geog, lon3, lat3)
        x_tm35fin, y_tm35fin, _, _ = etrs_exact_from_lonlat(triangles, vertices, lon3, lat3)
        assert abs(x_tm35fin_from_geog_interp - x_tm35fin) < 1e-3, (id3_geog, x_tm35fin_from_geog_interp, x_tm35fin)
        assert abs(y_tm35fin_from_geog_interp - y_tm35fin) < 1e-3, (id3_geog, y_tm35fin_from_geog_interp, y_tm35fin)

    is_in_shape = shape.Intersects(ogr.CreateGeometryFromWkt('POINT(%.18g %.18g)' % (lon1, lat1))) or \
                  shape.Intersects(ogr.CreateGeometryFromWkt('POINT(%.18g %.18g)' % (lon2, lat2))) or \
                  shape.Intersects(ogr.CreateGeometryFromWkt('POINT(%.18g %.18g)' % (lon3, lat3)))
    refinement_threshold = refinement_threshold_in_shape if is_in_shape else refinement_threshold_outside

    for w1, w2, w3 in ((1, 1, 1), (2, 1, 1), (1, 2, 1), (1, 1, 2), (4, 4, 1), (4, 1, 4), (1, 4, 4)):
        lon_middle = (w1 * lon1 + w2 * lon2 + w3 * lon3) / (w1 + w2 + w3)
        lat_middle = (w1 * lat1 + w2 * lat2 + w3 * lat3) / (w1 + w2 + w3)
        x_tm35fin_from_geog_interp, y_tm35fin_from_geog_interp, _, _ = etrs_approx_from_lonlat(vertices_geog, lon1, lat1, lon2, lat2, lon3, lat3, id1_geog, id2_geog, id3_geog, lon_middle, lat_middle)
        x_tm35fin, y_tm35fin, _, _ = etrs_exact_from_lonlat(triangles, vertices, lon_middle, lat_middle)
        square_error = (x_tm35fin_from_geog_interp - x_tm35fin)**2 + (y_tm35fin_from_geog_interp - y_tm35fin)**2
        error = square_error ** 0.5
        if error > refinement_threshold:
            break

    if error > refinement_threshold:
        print('Depth %d: Refine triangle (%s, %s, %s) due to error %f' % (depth, id1_geog, id2_geog, id3_geog, error))

        if False:
            # Refine by adding a point in centre and creating 3 triangles ==> does not converge
            id_middle = 'centre_of_' + id1_geog + '_' + id2_geog + '_' + id3_geog
            etrs_lon_exact, etrs_lat_exact, _ = etrs35fin_to_etrsgeog.TransformPoint(x_tm35fin, y_tm35fin)
            assert id_middle not in vertices_geog
            vertices_geog[id_middle] = [ lon_middle, lat_middle, etrs_lon_exact, etrs_lat_exact ]
            refine_triangle(vertices, id1_geog, id2_geog, id_middle, refinement_threshold, vertices_geog, new_triangles, map_long_to_short_id, depth+1)
            refine_triangle(vertices, id1_geog, id3_geog, id_middle, refinement_threshold, vertices_geog, new_triangles, map_long_to_short_id, depth+1)
            refine_triangle(vertices, id2_geog, id3_geog, id_middle, refinement_threshold, vertices_geog, new_triangles, map_long_to_short_id, depth+1)
            return

        x1, y1, _ = kkjgeog_to_kkjproj.TransformPoint(lon1, lat1)
        x2, y2, _ = kkjgeog_to_kkjproj.TransformPoint(lon2, lat2)
        x3, y3, _ = kkjgeog_to_kkjproj.TransformPoint(lon3, lat3)

        if False:
            # Refine by adding a point in the middle of each edge and creating 4 triangles.
            id_middle_12 = 'middle_of_' + min(id1_geog, id2_geog) + '_' + max(id1_geog, id2_geog)
            id_middle_13 = 'middle_of_' + min(id1_geog, id3_geog) + '_' + max(id1_geog, id3_geog)
            id_middle_23 = 'middle_of_' + min(id2_geog, id3_geog) + '_' + max(id2_geog, id3_geog)
            lon_middle_12 = (lon1 + lon2) / 2
            lat_middle_12 = (lat1 + lat2) / 2
            lon_middle_13 = (lon1 + lon3) / 2
            lat_middle_13 = (lat1 + lat3) / 2
            lon_middle_23 = (lon2 + lon3) / 2
            lat_middle_23 = (lat2 + lat3) / 2

            x_middle_12 = (x1 + x2) / 2
            y_middle_12 = (y1 + y2) / 2
            x_middle_13 = (x1 + x3) / 2
            y_middle_13 = (y1 + y3) / 2
            x_middle_23 = (x2 + x3) / 2
            y_middle_23 = (y2 + y3) / 2

            for lon, lat, x, y, id in [
                (lon_middle_12, lat_middle_12, x_middle_12, y_middle_12, id_middle_12),
                (lon_middle_13, lat_middle_13, x_middle_13, y_middle_13, id_middle_13),
                (lon_middle_23, lat_middle_23, x_middle_23, y_middle_23, id_middle_23) ]:
                if id not in vertices_geog:
                    try:
                        _, _, etrs_lon_exact, etrs_lat_exact = etrs_exact_from_lonlat(triangles, vertices, lon, lat)
                    except:
                        lon, lat, _ = kkjproj_to_kkjgeog.TransformPoint(x, y)
                        _, _, etrs_lon_exact, etrs_lat_exact = etrs_exact_from_lonlat(triangles, vertices, lon, lat)
                    vertices_geog[id] = [ lon, lat, etrs_lon_exact, etrs_lat_exact ]

            refine_triangle(vertices, id1_geog, id_middle_12, id_middle_13, refinement_threshold, vertices_geog, new_triangles, map_long_to_short_id, depth+1)
            refine_triangle(vertices, id2_geog, id_middle_12, id_middle_23, refinement_threshold, vertices_geog, new_triangles, map_long_to_short_id, depth+1)
            refine_triangle(vertices, id3_geog, id_middle_13, id_middle_23, refinement_threshold, vertices_geog, new_triangles, map_long_to_short_id, depth+1)
            refine_triangle(vertices, id_middle_12, id_middle_13, id_middle_23, refinement_threshold, vertices_geog, new_triangles)
            return

        # Refine by adding an intermediate point in the largest edge
        d12 = (x1 - x2)**2 + (y1 - y2) **2
        d13 = (x1 - x3)**2 + (y1 - y3) **2
        d23 = (x2 - x3)**2 + (y2 - y3) **2
        if d12 >= d13 and d12 >= d23:
            lon_middle = (lon1 + lon2) / 2
            lat_middle = (lat1 + lat2) / 2
            x_middle = (x1 + x2) / 2
            y_middle = (y1 + y2) / 2
            long_id = 'middle_of_' + min(id1_geog, id2_geog) + '_' + max(id1_geog, id2_geog)
            id_A1 = id1_geog
            id_A2 = id3_geog
            id_B1 = id2_geog
            id_B2 = id3_geog
        elif d13 >= d12 and d13 >= d23:
            lon_middle = (lon1 + lon3) / 2
            lat_middle = (lat1 + lat3) / 2
            x_middle = (x1 + x3) / 2
            y_middle = (y1 + y3) / 2
            long_id = 'middle_of_' + min(id1_geog, id3_geog) + '_' + max(id1_geog, id3_geog)
            id_A1 = id1_geog
            id_A2 = id2_geog
            id_B1 = id3_geog
            id_B2 = id2_geog
        else:
            lon_middle = (lon2 + lon3) / 2
            lat_middle = (lat2 + lat3) / 2
            x_middle = (x2 + x3) / 2
            y_middle = (y2 + y3) / 2
            long_id = 'middle_of_' + min(id2_geog, id3_geog) + '_' + max(id2_geog, id3_geog)
            id_A1 = id2_geog
            id_A2 = id1_geog
            id_B1 = id3_geog
            id_B2 = id1_geog

        if long_id not in map_long_to_short_id:
            try:
                _, _, etrs_lon_exact, etrs_lat_exact = etrs_exact_from_lonlat(triangles, vertices, lon_middle, lat_middle)
            except:
                lon_middle, lat_middle, _ = kkjproj_to_kkjgeog.TransformPoint(x_middle, y_middle)
                _, _, etrs_lon_exact, etrs_lat_exact = etrs_exact_from_lonlat(triangles, vertices, lon_middle, lat_middle)

            short_id = '_gen' + str(len(map_long_to_short_id)-1)
            assert short_id not in vertices_geog
            map_long_to_short_id[long_id] = short_id

            vertices_geog[short_id] = [ lon_middle, lat_middle, etrs_lon_exact, etrs_lat_exact ]

        short_id = map_long_to_short_id[long_id]

        refine_triangle(vertices, short_id, id_A1, id_A2, refinement_threshold_in_shape, shape, refinement_threshold_outside, vertices_geog, new_triangles, map_long_to_short_id, depth+1)
        refine_triangle(vertices, short_id, id_B1, id_B2, refinement_threshold_in_shape, shape, refinement_threshold_outside, vertices_geog, new_triangles, map_long_to_short_id, depth+1)

    else:
        new_triangles.append((id1_geog, id2_geog, id3_geog))


def refine_triangulation(vertices_geog, triangles, refinement_threshold_in_shape, shape, refinement_threshold_outside):

    new_vertices_geog = {}
    for k in vertices_geog:
        new_vertices_geog[k] = vertices_geog[k]
    new_triangles = []
    map_long_to_short_id = {}
    for (id1, id2, id3) in triangles:
        refine_triangle(vertices, id1, id2, id3, refinement_threshold_in_shape, shape, refinement_threshold_outside, new_vertices_geog, new_triangles, map_long_to_short_id)
    return (new_vertices_geog, new_triangles)

grid_min_easting, grid_max_easting, grid_stepx, grid_min_northing, grid_max_northing, grid_stepy, grid_delta_eastings = read_grid(os.path.join(dir_data, 'grid-YKJ-ETRSTM35FIN-L-I-1km-Windows.asc'))
grid_min_easting2, grid_max_easting2, grid_stepx2, grid_min_northing2, grid_max_northing2, grid_stepy2, grid_delta_northings = read_grid(os.path.join(dir_data, 'grid-YKJ-ETRSTM35FIN-P-E-1km-Windows.asc'))
assert grid_min_easting == grid_min_easting2
assert grid_max_easting == grid_max_easting2
assert grid_min_northing == grid_min_northing2
assert grid_max_northing == grid_max_northing2
assert grid_stepx == grid_stepx2
assert grid_stepy == grid_stepy2
grid_width = 1 + int((grid_max_easting - grid_min_easting) / grid_stepx + 0.5)
grid_height = 1 + int((grid_max_northing - grid_min_northing) / grid_stepy + 0.5)

vertices, vertices_geog = read_vertices()
triangles, envelope = read_triangles()
xmin, xmax, ymin, ymax = envelope.GetEnvelope()

def dump_triangles_geog(filename, layername, triangles, vertices):
    ds = ogr.GetDriverByName('GPKG').CreateDataSource(filename)
    lyr = ds.CreateLayer(layername, geom_type = ogr.wkbPolygon, srs = kkjgeog)
    for id1, id2, id3 in triangles:
        lon1 = vertices_geog[id1][0]
        lat1 = vertices_geog[id1][1]
        lon2 = vertices_geog[id2][0]
        lat2 = vertices_geog[id2][1]
        lon3 = vertices_geog[id3][0]
        lat3 = vertices_geog[id3][1]
        f = ogr.Feature(lyr.GetLayerDefn())
        tri_geom = ogr.CreateGeometryFromWkt('POLYGON((%.18g %.18g,%.18g %.18g,%.18g %.18g,%.18g %.18g))' % (lon1, lat1, lon2, lat2, lon3, lat3, lon1, lat1))
        f.SetGeometry(tri_geom)
        lyr.CreateFeature(f)

ds = ogr.Open(os.path.join(dir_data, 'finland.gpkg'))
lyr = ds.GetLayer(0)
f = lyr.GetNextFeature()
shape = f.GetGeometryRef().Clone()

if use_refinement:
    dump_triangles_geog(os.path.join(dir_data, 'triangles_geog.gpkg'), 'triangles_geog', triangles, vertices_geog)
    refinement_threshold_in_shape = 0.001
    refinement_threshold_outside = 0.01

    vertices_geog, triangles_geog = refine_triangulation(vertices_geog, triangles, refinement_threshold_in_shape, shape.Buffer(0.1), refinement_threshold_outside)
    dump_triangles_geog(os.path.join(dir_data, 'triangles_geog_refined.gpkg'), 'triangles_geog_refined', triangles_geog, vertices_geog)

else:
    triangles_geog = triangles

class MyFloat:
    def __init__(self, f):
        self.f = f

j = {}
j['vertices'] = []
j['triangles'] = []
for k in vertices_geog:
    j['vertices'].append([ k, MyFloat(vertices_geog[k][0]), MyFloat(vertices_geog[k][1]), MyFloat(vertices_geog[k][2]), MyFloat(vertices_geog[k][3]) ])
for id1, id2, id3 in triangles_geog:
    j['triangles'].append([id1, id2, id3])

class DecimalEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, MyFloat):
            return float('%.9f' % obj.f)
        return json.JSONEncoder.default(self, obj)

open(os.path.join(dir_data, 'triangulation.json'), 'wt').write(DecimalEncoder().encode(j))

print('New number of triangles: ', len(triangles_geog), 'Old number of triangles: ', len(triangles))
print('New number of vertices: ', len(vertices_geog), 'Old number of vertices: ', len(vertices))


# Check influence of northings in interpolation quality
#ymax = (1 * ymin + 2 * ymax) / 3

class Stat:
    def __init__(self):
        self.count_points = 0
        self.total_square_error = 0
        self.max_error = 0
        self.total_abs_error = 0
        self.points_below_1mm = 0
        self.points_below_2mm = 0
        self.points_below_3mm = 0
        self.points_below_4mm = 0
        self.points_below_5mm = 0
        self.points_below_1cm = 0
        self.points_below_2cm = 0
        self.points_below_5cm = 0
        self.points_below_10cm = 0
        self.points_below_20cm = 0
        self.total_error_easting = 0
        self.total_error_northing = 0
        self.total_error_lon_lat_valid = False
        self.total_error_lon = 0
        self.total_error_lat = 0

    def add_point(self, x_ref, y_ref, x, y, lon, lat):
        self.count_points += 1
        square_error = (x_ref - x)**2 + (y_ref - y)**2
        self.total_square_error += square_error
        error = square_error ** 0.5
        if error > self.max_error:
            self.max_error = error
        self.total_abs_error += error
        if error < 0.001:
            self.points_below_1mm += 1
        if error < 0.002:
            self.points_below_2mm += 1
        if error < 0.003:
            self.points_below_3mm += 1
        if error < 0.004:
            self.points_below_4mm += 1
        if error < 0.005:
            self.points_below_5mm += 1
        if error < 0.01:
            self.points_below_1cm += 1
        if error < 0.02:
            self.points_below_2cm += 1
        if error < 0.05:
            self.points_below_5cm += 1
        if error < 0.10:
            self.points_below_10cm += 1
        if error < 0.20:
            self.points_below_20cm += 1

        self.total_error_easting += x - x_ref
        self.total_error_northing += y - y_ref

        if not (lat == 0 and lon == 0):
            lon_ref, lat_ref, _ = etrs35fin_to_etrsgeog.TransformPoint(x_ref, y_ref)
            self.total_error_lon_lat_valid = True
            self.total_error_lon += lon - lon_ref
            self.total_error_lat += lat - lat_ref

    def display(self):
        print('RMSE: %f m' % (self.total_square_error / self.count_points)**0.5)
        print('Mean absolute error: %f m' % (self.total_abs_error / self.count_points))
        print('Max error: %f m' % (self.max_error))
        print('Percentage of points with max 1mm error: %.02f %%' % (100.0 * self.points_below_1mm / self.count_points))
        print('Percentage of points with max 2mm error: %.02f %%' % (100.0 * self.points_below_2mm / self.count_points))
        print('Percentage of points with max 3mm error: %.02f %%' % (100.0 * self.points_below_3mm / self.count_points))
        print('Percentage of points with max 4mm error: %.02f %%' % (100.0 * self.points_below_4mm / self.count_points))
        print('Percentage of points with max 5mm error: %.02f %%' % (100.0 * self.points_below_5mm / self.count_points))
        print('Percentage of points with max 1cm error: %.02f %%' % (100.0 * self.points_below_2cm / self.count_points))
        print('Percentage of points with max 2cm error: %.02f %%' % (100.0 * self.points_below_2cm / self.count_points))
        print('Percentage of points with max 5cm error: %.02f %%' % (100.0 * self.points_below_5cm / self.count_points))
        print('Percentage of points with max 10cm error: %.02f %%' % (100.0 * self.points_below_10cm / self.count_points))
        print('Percentage of points with max 20cm error: %.02f %%' % (100.0 * self.points_below_20cm / self.count_points))

        bias_easting = self.total_error_easting / self.count_points
        bias_northing = self.total_error_northing / self.count_points
        print('Bias in eastings: %f m' % bias_easting)
        print('Bias in northings: %f m' % bias_northing)

        if self.total_error_lon_lat_valid:
            bias_lon = self.total_error_lon / self.count_points
            bias_lat = self.total_error_lat / self.count_points
            print('Bias in longitude: %f deg' % bias_lon)
            print('Bias in latitude: %f deg' % bias_lat)


stat_all = Stat()
stat_inside = Stat()
stat_outside = Stat()
stat_grid = Stat()

f = open('output.csv', 'wt')
f.write('easting_kkj,northing_kkj,easting_etrs_exact,northing_etrs_exact,easting_etrs_interpolated,northing_etrs_interpolated,easting_interpolated_minus_exact,northing_interpolated_minus_exact\n')

def process_point(x,y):

    kkj_lon, kkj_lat, _ = kkjproj_to_kkjgeog.TransformPoint(x, y)

    p = ogr.CreateGeometryFromWkt('POINT (%.18g %.18g)' % (x,y))
    if not p.Intersects(envelope):
        return

    # Perform barycentric interpolation on projected coordinates
    ret = find_triangle(triangles, vertices, x, y)
    assert ret
    id1, id2, id3, lambda_1, lambda_2, lambda_3 = ret
    x_tm35fin = vertices[id1][2] * lambda_1 + vertices[id2][2] * lambda_2 + vertices[id3][2] * lambda_3
    y_tm35fin = vertices[id1][3] * lambda_1 + vertices[id2][3] * lambda_2 + vertices[id3][3] * lambda_3

    if x >= grid_min_easting and x <= grid_max_easting and y >= grid_min_northing and y <= grid_max_northing:
        float_i = (x - grid_min_easting) / grid_stepx
        i = int(float_i)
        float_j = (y - grid_min_northing) / grid_stepy
        j = int(float_j)
        if i + 1 < grid_width and j + 1 < grid_height and \
            grid_delta_eastings[j][i] != -999999.0 and \
            grid_delta_eastings[j][i+1] != -999999.0 and \
            grid_delta_eastings[j+1][i] != -999999.0 and \
            grid_delta_eastings[j+1][i+1] != -999999.0 and \
            grid_delta_northings[j][i] != -999999.0 and \
            grid_delta_northings[j][i+1] != -999999.0 and \
            grid_delta_northings[j+1][i] != -999999.0 and \
            grid_delta_northings[j+1][i+1] != -999999.0:

            di = float_i - i
            dj = float_j - j
            assert di >= 0 and di <= 1
            assert dj >= 0 and dj <= 1
            dx = (1 - di) * (1 - dj) * grid_delta_eastings[j][i] + \
                    di * (1 - dj) * grid_delta_eastings[j][i+1] + \
                    (1 - di) * dj * grid_delta_eastings[j+1][i] + \
                    di * dj * grid_delta_eastings[j+1][i+1]
            dy = (1 - di) * (1 - dj) * grid_delta_northings[j][i] + \
                    di * (1 - dj) * grid_delta_northings[j][i+1] + \
                    (1 - di) * dj * grid_delta_northings[j+1][i] + \
                    di * dj * grid_delta_northings[j+1][i+1]

            x_tm35fin_grid = x + dx
            y_tm35fin_grid = y + dy
            # print (x_tm35fin, y_tm35fin, x_tm35fin_grid, y_tm35fin_grid)

            stat_grid.add_point(x_tm35fin, y_tm35fin, x_tm35fin_grid, y_tm35fin_grid, 0, 0)

    # Perform barycentric interpolation on geographic coordinates
    ret = find_triangle(triangles_geog, vertices_geog, kkj_lon, kkj_lat)
    if ret is None:
        print('Skipping ', (x, y, kkj_lon, kkj_lat))
        return
    id1, id2, id3, lambda_1, lambda_2, lambda_3 = ret
    etrs_lon = vertices_geog[id1][2] * lambda_1 + vertices_geog[id2][2] * lambda_2 + vertices_geog[id3][2] * lambda_3
    etrs_lat = vertices_geog[id1][3] * lambda_1 + vertices_geog[id2][3] * lambda_2 + vertices_geog[id3][3] * lambda_3
    x_tm35fin_from_geog_interp, y_tm35fin_from_geog_interp, _ = etrsgeog_to_etrs35fin.TransformPoint(etrs_lon, etrs_lat)

    #print(x, y, x_tm35fin, y_tm35fin, x_tm35fin_from_geog_interp, y_tm35fin_from_geog_interp)

    f.write('%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,%.15g\n' % (x, y, x_tm35fin, y_tm35fin, x_tm35fin_from_geog_interp, y_tm35fin_from_geog_interp, x_tm35fin_from_geog_interp - x_tm35fin, y_tm35fin_from_geog_interp - y_tm35fin))

    stat_all.add_point(x_tm35fin, y_tm35fin, x_tm35fin_from_geog_interp, y_tm35fin_from_geog_interp, etrs_lon, etrs_lat)

    if shape.Intersects(ogr.CreateGeometryFromWkt('POINT (%.18g %.18g)' % (kkj_lon,kkj_lat))):
        stat_inside.add_point(x_tm35fin, y_tm35fin, x_tm35fin_from_geog_interp, y_tm35fin_from_geog_interp, etrs_lon, etrs_lat)
    else:
        stat_outside.add_point(x_tm35fin, y_tm35fin, x_tm35fin_from_geog_interp, y_tm35fin_from_geog_interp, etrs_lon, etrs_lat)

# Sample random points inside the triangulation
#for i in range(20000):
#    x = random.uniform(xmin, xmax)
#    y = random.uniform(ymin, ymax)
#    process_point(x, y)

#for k in vertices:
#    x = vertices[k][0]
#    y = vertices[k][1]
#    process_point(x, y)

for (id1, id2, id3) in triangles:
    x1 = vertices[id1][0]
    y1 = vertices[id1][1]
    x2 = vertices[id2][0]
    y2 = vertices[id2][1]
    x3 = vertices[id3][0]
    y3 = vertices[id3][1]
    i = 0
    while i < 10:
        lambda_1 = random.uniform(0, 1)
        lambda_2 = random.uniform(0, 1)
        if lambda_1 + lambda_2 > 1:
            continue
        i += 1
        lambda_3 = 1 - lambda_1 - lambda_2
        x = lambda_1 * x1 + lambda_2 * x2 + lambda_3 * x3
        y = lambda_1 * y1 + lambda_2 * y2 + lambda_3 * y3
        process_point(x, y)

print('-------------------------------------------------------------------------------------------------------------')
print('Triangulation from geographic coordinates interpolation compared to triangulation from projected coordinates:')
print('-------------------------------------------------------------------------------------------------------------')
print('\n')
print('- Global stats:')
stat_all.display()
print('\n')
print('- Points inside shape:')
stat_inside.display()
print('\n')
print('- Points outside shape:')
stat_outside.display()

print('\n')
print('------------------------------------------------------------------------')
print('Grid interpolation compared to triangulation from projected coordinates:')
print('------------------------------------------------------------------------')
stat_grid.display()

f.close()
