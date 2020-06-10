x1 = 1.0
y1 = -2.0
z1 = 0.0

x2 = 3.0
y2 = 1.0
z2 = 4.0

x3 = 0.0
y3 = -1.0
z3 = 2.0

import numpy

x = numpy.array([[x1, y1, 1],[x2, y2, 1],[x3, y3, 1]])
#M = numpy.matmul(numpy.linalg.inv(x), numpy.array([[z1],[z2],[z3]]))
M = numpy.linalg.solve(x, numpy.array([[z1],[z2],[z3]]))  # slightly more exact

def plane_z(x, y):
    return M[0][0] * x + M[1][0] * y + M[2][0]

print(plane_z(x1,y1))
print(plane_z(x2,y2))
print(plane_z(x3,y3))
print(plane_z(0,0))
print(plane_z((x1 + 2*x2) /3, (y1+2*y2)/3))
print(plane_z((x1 + x2 + x3)/3, (y1+y2+y3)/3))

print('')

def barycentric_interp(x, y):
    det_T = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)
    lambda_1 = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / det_T
    lambda_2 = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / det_T
    lambda_3 = 1 - lambda_1 - lambda_2
    return lambda_1 * z1 + lambda_2 * z2 + lambda_3 * z3


print(barycentric_interp(x1,y1))
print(barycentric_interp(x2,y2))
print(barycentric_interp(x3,y3))
print(barycentric_interp(0,0))
print(barycentric_interp((x1 + 2*x2) /3, (y1+2*y2)/3))
print(barycentric_interp((x1 + x2 + x3)/3, (y1+y2+y3)/3))