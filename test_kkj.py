from osgeo import ogr, osr
import random
import os

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
#kkjgeog_to_kkjproj = osr.CoordinateTransformation(kkjgeog, kkjproj)

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
        if lambda_1 >= 0 and lambda_1 <= 1 and \
           lambda_2 >= 0 and lambda_2 <= 1 and \
           lambda_3 >= 0 and lambda_3 <= 1:
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

# Check influence of northings in interpolation quality
#ymax = (1 * ymin + 2 * ymax) / 3

bias_easting = 0
bias_northing = 0

for iteration in range(2):

    if iteration == 1:
        print('\nDoing a second iteration with bias compensation:')

    count_points = 0
    total_square_error = 0
    max_error = 0
    total_abs_error = 0
    points_below_2cm = 0
    points_below_5cm = 0
    points_below_10cm = 0
    points_below_20cm = 0
    total_error_easting = 0
    total_error_northing = 0

    count_points_grid = 0
    total_square_error_grid = 0
    total_abs_error_grid = 0
    max_error_grid = 0
    total_error_easting_grid = 0
    total_error_northing_grid = 0

    # Sample random points inside the triangulation
    for i in range(20000):
        x = random.uniform(xmin, xmax)
        y = random.uniform(ymin, ymax)

    #for k in vertices:
    #    x = vertices[k][0]
    #    y = vertices[k][1]
        p = ogr.CreateGeometryFromWkt('POINT (%.18g %.18g)' % (x,y))
        if not p.Intersects(envelope):
            continue

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
                count_points_grid += 1
                square_error = (x_tm35fin_grid - x_tm35fin)**2 + (y_tm35fin_grid - y_tm35fin)**2
                total_square_error_grid += square_error
                error = square_error ** 0.5
                if error > max_error_grid:
                    max_error_grid = error
                total_abs_error_grid += error

                total_error_easting_grid += x_tm35fin_grid - x_tm35fin
                total_error_northing_grid += y_tm35fin_grid - y_tm35fin

        # Perform barycentric interpolation on geographic coordinates
        kkj_lon, kkj_lat, _ = kkjproj_to_kkjgeog.TransformPoint(x, y)
        ret = find_triangle(triangles, vertices_geog, kkj_lon, kkj_lat)
        if ret is None:
            print('Skipping ', (x, y, kkj_lon, kkj_lat))
            continue
        id1, id2, id3, lambda_1, lambda_2, lambda_3 = ret
        etrs_lon = vertices_geog[id1][2] * lambda_1 + vertices_geog[id2][2] * lambda_2 + vertices_geog[id3][2] * lambda_3
        etrs_lat = vertices_geog[id1][3] * lambda_1 + vertices_geog[id2][3] * lambda_2 + vertices_geog[id3][3] * lambda_3
        x_tm35fin_from_geog_interp, y_tm35fin_from_geog_interp, _ = etrsgeog_to_etrs35fin.TransformPoint(etrs_lon, etrs_lat)

        x_tm35fin_from_geog_interp -= bias_easting
        y_tm35fin_from_geog_interp -= bias_northing

        #print(x, y, x_tm35fin, y_tm35fin, x_tm35fin_from_geog_interp, y_tm35fin_from_geog_interp)
        count_points += 1
        square_error = (x_tm35fin_from_geog_interp - x_tm35fin)**2 + (y_tm35fin_from_geog_interp - y_tm35fin)**2
        total_square_error += square_error
        error = square_error ** 0.5
        if error > max_error:
            max_error = error
        total_abs_error += error
        if error < 0.02:
            points_below_2cm += 1
        if error < 0.05:
            points_below_5cm += 1
        if error < 0.10:
            points_below_10cm += 1
        if error < 0.20:
            points_below_20cm += 1

        total_error_easting += x_tm35fin_from_geog_interp - x_tm35fin
        total_error_northing += y_tm35fin_from_geog_interp - y_tm35fin

    print('Triangulation from geographic coordinates interpolation compared to triangulation from projected coordinates:')
    print('RMSE: %f m' % (total_square_error / count_points)**0.5)
    print('Mean absolute error: %f m' % (total_abs_error / count_points))
    print('Max error: %f m' % (max_error))
    print('Percentage of points with max 2cm error: %.02f %%' % (100.0 * points_below_2cm / count_points))
    print('Percentage of points with max 5cm error: %.02f %%' % (100.0 * points_below_5cm / count_points))
    print('Percentage of points with max 10cm error: %.02f %%' % (100.0 * points_below_10cm / count_points))
    print('Percentage of points with max 20cm error: %.02f %%' % (100.0 * points_below_20cm / count_points))

    bias_easting = total_error_easting / count_points
    bias_northing = total_error_northing / count_points

    print('Bias in eastings: %f' % bias_easting)
    print('Bias in northings: %f' % bias_northing)

    if iteration == 0:
        print('\n')
        print('Grid interpolation compared to triangulation from projected coordinates:')
        print('RMSE: %f m' % (total_square_error_grid / count_points_grid)**0.5)
        print('Mean absolute error: %f m' % (total_abs_error_grid / count_points_grid))
        print('Max error: %f m' % (max_error_grid))
        print('Bias in eastings: %f' % (total_error_easting_grid / count_points_grid))
        print('Bias in northings: %f' % (total_error_northing_grid / count_points_grid))
