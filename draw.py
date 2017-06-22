from display import *
from matrix import *
from math import *
from gmath import *
from sys import maxint

def add_polygon( polygons, x0, y0, z0, x1, y1, z1, x2, y2, z2 ):
    add_point(polygons, x0, y0, z0);
    add_point(polygons, x1, y1, z1);
    add_point(polygons, x2, y2, z2);

def draw_polygons( matrix, screen, zbuffer, color ):
    if len(matrix) < 2:
        print 'Need at least 3 points to draw'
        return 

    point = 0    
    while point < len(matrix) - 2:

        normal = calculate_normal(matrix, point)[:]
        #print normal


        if normal[2] > 0:
            scanline(matrix, point, screen, color, zbuffer)
            draw_line( int(matrix[point][0]),
                       int(matrix[point][1]),
                       
                       zbuffer[int(matrix[point][0])][int(matrix[point][1])], 
                       #matrix[point][2],
                       int(matrix[point+1][0]),
                       int(matrix[point+1][1]),
                       zbuffer[int(matrix[point+1][0])][int(matrix[point+1][1])],                       
                       #matrix[point+1][2],
                       screen, zbuffer, color)
            print [zbuffer[int(matrix[point][0])][int(matrix[point][1])], matrix[point][2]]
            draw_line( int(matrix[point+2][0]),
                       int(matrix[point+2][1]),
                       zbuffer[int(matrix[point+2][0])][int(matrix[point+2][1])], 
                       #matrix[point+2][2],
                       int(matrix[point+1][0]),
                       int(matrix[point+1][1]),
                       zbuffer[int(matrix[point+1][0])][int(matrix[point+1][1])], 
                       #matrix[point+1][2],
                       screen, zbuffer, color)
            draw_line( int(matrix[point][0]),
                       int(matrix[point][1]),
                       zbuffer[int(matrix[point][0])][int(matrix[point][1])], 
                       #matrix[point][2],
                       int(matrix[point+2][0]),
                       int(matrix[point+2][1]),
                       zbuffer[int(matrix[point+2][0])][int(matrix[point+2][1])], 
                       #matrix[point+2][2],
                       screen, zbuffer, color)    
        point+= 3

#given the bottom, mid and top values of y, get the matching x for each
def findMatchingx(matrix, point, By, My, Ty):
    Botx = [maxint, maxint, maxint]
    Botz = [maxint, maxint, maxint]
    Topx = [-maxint - 1, -maxint - 1, -maxint - 1]
    Topz = [-maxint - 1, -maxint - 1, -maxint - 1]
    ans = [0, 0, 0, 0, 0, 0] #ans = [Botx, Midx, Topx, Botz, Midz, Topz]
    #collect all possible values of Bx and Tx, Bz and Tz
    for i in range(3):
        if matrix[point+i][1] == By:
            Botx[i] = matrix[point+i][0]
            Botz[i] = matrix[point+i][2]
    for i in range(3):
        if matrix[point+i][1] == Ty:
            Topx[i] = matrix[point+i][0]
            Topz[i] = matrix[point+i][2]
    ans[0] = min(Botx[0], Botx[1], Botx[2])
    ans[2] = max(Topx[0], Topx[1], Topx[2])
    ans[1] = (matrix[point][0] + matrix[point+1][0] + matrix[point+2][0]) - (ans[0] + ans[2])
    ans[3] = min(Botz[0], Botz[1], Botz[2])
    ans[5] = max(Topz[0], Topz[1], Topz[2])
    ans[4] = (matrix[point][2] + matrix[point+1][2] + matrix[point+2][2]) - (ans[3] + ans[5])
    return ans

def scanline(matrix, point, screen, color, zbuffer):
    color = [0, 255, 255]
    increment_constant = 1
    x1 = matrix[point][0]
    y1 = matrix[point][1]
    x2 = matrix[point+1][0]
    y2 = matrix[point+1][1]
    x3 = matrix[point+2][0]
    y3 = matrix[point+2][1]
    By = min(y1, y2, y3)
    Ty = max(y1, y2, y3)
    My = y1 + y2 + y3 - (By + Ty)
    xzValues = findMatchingx(matrix, point, By, My, Ty)
    Bx = xzValues[0]
    Mx = xzValues[1]
    Tx = xzValues[2]
    Bz = xzValues[3]
    Mz = xzValues[4]
    Tz = xzValues[5]
    if Ty-By == 0:
        slopeTopx = 0
        slopeTopz = 0
    else:
        slopeTopx = (Tx - Bx) * 1.0 / (Ty - By)
        slopeTopz = (Tz - Bz) * 1.0 / (Ty - By)
    if My-By == 0:
        slopeMidx = 0
        slopeMidz = 0
    else:
        slopeMidx = (Mx - Bx) * 1.0 / (My - By)
        slopeMidz = (Mz - Bz) * 1.0 / (My - By)
    y = By
    xTop = Bx
    xMid = Bx
    zTop = Bz
    zMid = Bz
    #split the triangle to two parts with horizontal bases
    
    while y <= My-.3:
        xTop = xTop + slopeTopx * increment_constant
        xMid = xMid + slopeMidx * increment_constant
        zTop = zTop + slopeTopz * increment_constant
        zMid = zMid + slopeMidz * increment_constant
        draw_line(int(xMid), int(y), int(zMid),
                  int(xTop), int(y), int(zTop),
                  screen, zbuffer, color)
        y += increment_constant

    xMid = Mx
    if Ty-My == 0:
        slopeMidx = 0
        slopeMidz = 0
    else:
        slopeMidx = (Tx - Mx) * 1.0 / (Ty - My)
        slopeMidz = (Tz - Mz) * 1.0 / (Ty - My)
    while y <= Ty-.3:
        xTop = xTop + slopeTopx * increment_constant
        xMid = xMid + slopeMidx * increment_constant
        zTop = zTop + slopeTopz * increment_constant
        zMid = zMid + slopeMidz * increment_constant
        draw_line(int(xTop), int(y), int(zTop),
                  int(xMid), int(y), int(zMid),
                 screen, zbuffer, color)
        y += increment_constant


def add_box( polygons, x, y, z, width, height, depth ):
    x1 = x + width
    y1 = y - height
    z1 = z - depth

    #front
    add_polygon(polygons, x, y, z, x1, y1, z, x1, y, z);
    add_polygon(polygons, x, y, z, x, y1, z, x1, y1, z);
  
    #back
    add_polygon(polygons, x1, y, z1, x, y1, z1, x, y, z1);
    add_polygon(polygons, x1, y, z1, x1, y1, z1, x, y1, z1);
  
    #right side
    add_polygon(polygons, x1, y, z, x1, y1, z1, x1, y, z1);
    add_polygon(polygons, x1, y, z, x1, y1, z, x1, y1, z1);
    #left side
    add_polygon(polygons, x, y, z1, x, y1, z, x, y, z);
    add_polygon(polygons, x, y, z1, x, y1, z1, x, y1, z);
  
    #top
    add_polygon(polygons, x, y, z1, x1, y, z, x1, y, z1);
    add_polygon(polygons, x, y, z1, x, y, z, x1, y, z);
    #bottom
    add_polygon(polygons, x, y1, z, x1, y1, z1, x1, y1, z);
    add_polygon(polygons, x, y1, z, x, y1, z1, x1, y1, z1);

def add_sphere( edges, cx, cy, cz, r, step ):
    points = generate_sphere(cx, cy, cz, r, step)
    num_steps = int(1/step+0.1)
    
    lat_start = 0
    lat_stop = num_steps
    longt_start = 0
    longt_stop = num_steps

    num_steps+= 1
    for lat in range(lat_start, lat_stop):
        for longt in range(longt_start, longt_stop):
            
            p0 = lat * (num_steps) + longt
            p1 = p0+1
            p2 = (p1+num_steps) % (num_steps * (num_steps-1))
            p3 = (p0+num_steps) % (num_steps * (num_steps-1))

            if longt != num_steps - 2:
	        add_polygon( edges, points[p0][0],
		             points[p0][1],
		             points[p0][2],
		             points[p1][0],
		             points[p1][1],
		             points[p1][2],
		             points[p2][0],
		             points[p2][1],
		             points[p2][2])
            if longt != 0:
	        add_polygon( edges, points[p0][0],
		             points[p0][1],
		             points[p0][2],
		             points[p2][0],
		             points[p2][1],
		             points[p2][2],
		             points[p3][0],
		             points[p3][1],
		             points[p3][2])

def generate_sphere( cx, cy, cz, r, step ):
    points = []
    num_steps = int(1/step+0.1)
    
    rot_start = 0
    rot_stop = num_steps
    circ_start = 0
    circ_stop = num_steps
            
    for rotation in range(rot_start, rot_stop):
        rot = step * rotation
        for circle in range(circ_start, circ_stop+1):
            circ = step * circle

            x = r * math.cos(math.pi * circ) + cx
            y = r * math.sin(math.pi * circ) * math.cos(2*math.pi * rot) + cy
            z = r * math.sin(math.pi * circ) * math.sin(2*math.pi * rot) + cz

            points.append([x, y, z])
            #print 'rotation: %d\tcircle%d'%(rotation, circle)
    return points
        
def add_torus( edges, cx, cy, cz, r0, r1, step ):
    points = generate_torus(cx, cy, cz, r0, r1, step)
    num_steps = int(1/step+0.1)
    
    lat_start = 0
    lat_stop = num_steps
    longt_start = 0
    longt_stop = num_steps
    
    for lat in range(lat_start, lat_stop):
        for longt in range(longt_start, longt_stop):

            p0 = lat * (num_steps) + longt;
            if (longt == num_steps - 1):
	        p1 = p0 - longt;
            else:
	        p1 = p0 + 1;
            p2 = (p1 + num_steps) % (num_steps * num_steps);
            p3 = (p0 + num_steps) % (num_steps * num_steps);

            add_polygon(edges,
                        points[p0][0],
                        points[p0][1],
                        points[p0][2],
                        points[p3][0],
                        points[p3][1],
                        points[p3][2],
                        points[p2][0],
                        points[p2][1],
                        points[p2][2] )
            add_polygon(edges,
                        points[p0][0],
                        points[p0][1],
                        points[p0][2],
                        points[p2][0],
                        points[p2][1],
                        points[p2][2],
                        points[p1][0],
                        points[p1][1],
                        points[p1][2] )

def generate_torus( cx, cy, cz, r0, r1, step ):
    points = []
    num_steps = int(1/step+0.1)
    
    rot_start = 0
    rot_stop = num_steps
    circ_start = 0
    circ_stop = num_steps
    
    for rotation in range(rot_start, rot_stop):
        rot = step * rotation
        for circle in range(circ_start, circ_stop):
            circ = step * circle

            x = math.cos(2*math.pi * rot) * (r0 * math.cos(2*math.pi * circ) + r1) + cx;
            y = r0 * math.sin(2*math.pi * circ) + cy;
            z = -1*math.sin(2*math.pi * rot) * (r0 * math.cos(2*math.pi * circ) + r1) + cz;

            points.append([x, y, z])
    return points

def add_circle( points, cx, cy, cz, r, step ):
    x0 = r + cx
    y0 = cy
    t = step

    while t <= 1.00001:
        x1 = r * math.cos(2*math.pi * t) + cx;
        y1 = r * math.sin(2*math.pi * t) + cy;

        add_edge(points, x0, y0, cz, x1, y1, cz)
        x0 = x1
        y0 = y1
        t+= step

def add_curve( points, x0, y0, x1, y1, x2, y2, x3, y3, step, curve_type ):

    xcoefs = generate_curve_coefs(x0, x1, x2, x3, curve_type)[0]
    ycoefs = generate_curve_coefs(y0, y1, y2, y3, curve_type)[0]

    t = step
    while t <= 1.00001:
        x = xcoefs[0] * t*t*t + xcoefs[1] * t*t + xcoefs[2] * t + xcoefs[3]
        y = ycoefs[0] * t*t*t + ycoefs[1] * t*t + ycoefs[2] * t + ycoefs[3]
                
        add_edge(points, x0, y0, 0, x, y, 0)
        x0 = x
        y0 = y
        t+= step

def draw_lines( matrix, screen, zbuffer, color ):
    if len(matrix) < 2:
        print 'Need at least 2 points to draw'
        return
    
    point = 0
    while point < len(matrix) - 1:
        draw_line( int(matrix[point][0]),
                   int(matrix[point][1]),
                   matrix[point][2],
                   int(matrix[point+1][0]),
                   int(matrix[point+1][1]),
                   matrix[point+1][2],
                   screen, zbuffer, color)    
        point+= 2
        
def add_edge( matrix, x0, y0, z0, x1, y1, z1 ):
    add_point(matrix, x0, y0, z0)
    add_point(matrix, x1, y1, z1)
    
def add_point( matrix, x, y, z=0 ):
    matrix.append( [x, y, z, 1] )
    



def draw_line( x0, y0, z0, x1, y1, z1, screen, zbuffer, color ):

    #swap points if going right -> left
    if x0 > x1:
        xt = x0
        yt = y0
        zt = z0
        x0 = x1
        y0 = y1
        z0 = z1
        x1 = xt
        y1 = yt
        z1 = zt

    x = x0
    y = y0
    z = z0
    A = 2 * (y1 - y0)
    B = -2 * (x1 - x0)
    wide = False
    tall = False

    if ( abs(x1-x0) >= abs(y1 - y0) ): #octants 1/8
        wide = True
        loop_start = x
        loop_end = x1
        dx_east = dx_northeast = 1
        dy_east = 0
        d_east = A
        distance = x1 - x
        if ( A > 0 ): #octant 1
            d = A + B/2
            dy_northeast = 1
            d_northeast = A + B
        else: #octant 8
            d = A - B/2
            dy_northeast = -1
            d_northeast = A - B

    else: #octants 2/7
        tall = True
        dx_east = 0
        dx_northeast = 1
        distance = abs(y1 - y)
        if ( A > 0 ): #octant 2
            d = A/2 + B
            dy_east = dy_northeast = 1
            d_northeast = A + B
            d_east = B
            loop_start = y
            loop_end = y1
        else: #octant 7
            d = A/2 - B
            dy_east = dy_northeast = -1
            d_northeast = A - B
            d_east = -1 * B
            loop_start = y1
            loop_end = y
    if (loop_start - loop_end != 0):
        dz = float(max(z0, z1) - min(z0, z1)) / (loop_start - loop_end)
    else:
        dz = 0

    while ( loop_start < loop_end ):
        #print [z0, z, z1]
        if z > zbuffer[x][y]:
            plot( screen, zbuffer, color, x, y, z )
        if ( (wide and ((A > 0 and d > 0) or (A < 0 and d < 0))) or
             (tall and ((A > 0 and d < 0) or (A < 0 and d > 0 )))):
            x+= dx_northeast
            y+= dy_northeast
            d+= d_northeast
        else:
            x+= dx_east
            y+= dy_east
            d+= d_east
        z += dz
        loop_start+= 1

    plot( screen, zbuffer, color, x, y, z )

    
