import numpy

# Triangle 1278
#1	6718527.414	3106266.213	6715706.377	106256.36
#391	6732415.306	3127183.015	6729588.768	127164.654
#669	6764362.883	3064653.319	6761522.932	64660.087

a1=-0.0000041547146677
a2=0.9995960042236450
deltaE=-2998727.0210
b1=0.9995935741323370
b2=0.0000068555931313
deltaN=-111.7490

x1 = 3106266.213
y1 = 6718527.414
x2 = 3127183.015
y2 = 6732415.306
x3 = 3064653.319
y3 = 6764362.883

X1 = 106256.36
Y1 = 6715706.377
X2 = 127164.654
Y2 = 6729588.768
X3 = 64660.087
Y3 = 6761522.932

# Recompute a1, a2, etc...
x = numpy.array([[x1, y1, 1],[x2, y2, 1],[x3, y3, 1]])
M = numpy.matmul(numpy.linalg.inv(x), numpy.array([[X1],[X2],[X3]]))
#M = numpy.linalg.solve(x, numpy.array([[X1],[X2],[X3]]))  # slightly more exact
a2 = M[0][0]
a1 = M[1][0]
deltaE = M[2][0]
M = numpy.matmul(numpy.linalg.inv(x), numpy.array([[Y1],[Y2],[Y3]]))
#M = numpy.linalg.solve(x, numpy.array([[Y1],[Y2],[Y3]]))  # slightly more exact
b2 = M[0][0]
b1 = M[1][0]
deltaN = M[2][0]
print(a1,a2, deltaE, b1, b2, deltaN)

def linear_2D(x,y):
    return deltaE + a1 * y + a2 * x, deltaN + b1 * y + b2 * x

def get_barycentric_coord(x1, y1, x2, y2, x3, y3, x, y):
    det_T = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)
    lambda_1 = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / det_T
    lambda_2 = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / det_T
    lambda_3 = 1 - lambda_1 - lambda_2
    return lambda_1, lambda_2, lambda_3

def interp_barycentric(x,y):
    lambda_1, lambda_2, lambda_3 = get_barycentric_coord(x1, y1, x2, y2, x3, y3, x, y)
    return lambda_1 * X1 + lambda_2 * X2 + lambda_3 * X3, \
           lambda_1 * Y1 + lambda_2 * Y2 + lambda_3 * Y3

print('With linear 2D interpolation:')
print(linear_2D(x1, y1))
print(linear_2D(x2, y2))
print(linear_2D(x3, y3))
print(linear_2D((x1+x2+x3)/3, (y1+y2+y3)/3))
print(linear_2D((x1+x2+2*x3)/4, (y1+y2+2*y3)/4))
print(linear_2D(0,0))

print('With barycentric interpolation:')
print(interp_barycentric(x1, y1))
print(interp_barycentric(x2, y2))
print(interp_barycentric(x3, y3))
print(interp_barycentric((x1+x2+x3)/3, (y1+y2+y3)/3))
print(interp_barycentric((x1+x2+2*x3)/4, (y1+y2+2*y3)/4))
print(interp_barycentric(0,0))