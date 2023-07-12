"""
This program simulates diffusion limited aggregation on an LxL grid.
Particles are initiated until the centre point is filled.
Author: Nico Grisouard, University of Toronto
Based on Paul J Kushner's DAL-eample.py
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc


def nextmove(x, y):
    """ randomly choose a direction
    1 = up, 2 = down, 3 = left, 4 = right"""
    direction =  np.random.randint(1, 5)
    if direction == 1:  # move up
        y += 1
    elif direction == 2:  # move down
        y -= 1
    elif direction == 3:  # move right
        x += 1
    elif direction == 4:  # move left
        x -= 1
    else:
        print("error: direction isn't 1-4")

    return x, y

def neighbour(x, y):
    """Return coords immidiately adj. to the cell
        |X|
       X|A|X for any cell A
        |X|
    """
    x_right = min(x+1, Lp-1)
    x_left = max(x-1, 0)
    y_up = min(y+1, Lp-1)
    y_down = max(y-1, 0)
    return [[x_right, y],[x_left, y],[x, y_up], [x, y_down]]

# %% main program starts here ------------------------------------------------|
# YOU NEED TO FINISH IT!


Lp = 101  # size of domain
N = 0  # number of particles (start count on 0)
# array to represent whether each gridpoint has an anchored particle
anchored = np.zeros((Lp, Lp), dtype=int)
# list to represent x and y positions of anchored points
anchored_points = [[], []]

centre_point = (Lp-1)//2  # middle point of domain

#Keep releasing particle unless dead
#Added a 'cap' of 1000 to save time
dead = False
while dead == False:    
    xp = centre_point
    yp = centre_point
    i = 0  # counter to keep track of animation of moving particle
    match = False
    
    #do while loop to update particle position until anchored
    while xp in range(1, Lp-1) and yp in range(1, Lp-1) and match == False:
        #new move for each step
        xpp, ypp = nextmove(xp, yp)
        #round position on grid to prevent index error
        xpp = max(min(Lp, xpp), 0)
        ypp = max(min(Lp, ypp), 0)
        #check if perspective position contacts with other anchored cells
        cells_adj = neighbour(xpp, ypp)
        for cell in cells_adj:
            if anchored[cell[1]][cell[0]] == 1:
                match = True
        xp = np.copy(xpp)
        yp = np.copy(ypp)
    #anchor the final position to the wall (or other particles)
    anchored_points[0].append(xp)
    anchored_points[1].append(yp)
    anchored[int(yp)][int(xp)] = 1
    N += 1
    print(N)
    #check if center point is anchored.
    if anchored[centre_point][centre_point] == 1:
        dead = True
    elif N == 1000:
        dead = True
    if dead == True:
        break
        
# set up animation of anchored points
#plot the walls
bottom = np.zeros(Lp)
top = bottom + Lp-1
x_box = np.arange(0, Lp)
left = np.zeros(Lp)
right = left + Lp-1
y_box = np.arange(0, Lp)
#plot the particles
plt.figure(figsize = (6, 4))
plt.title('DLA run for {} particles'.format(N))
plt.plot(anchored_points[0], anchored_points[1], '.', markersize=5)
plt.plot(x_box, bottom, color = 'blue')
plt.plot(x_box, top, color = 'blue')
plt.plot(top, y_box, color = 'blue')
plt.plot(bottom, y_box, color = 'blue')
plt.xlim([-10, Lp+10])
plt.ylim([-10, Lp+10])
plt.axis('equal')
plt.xlabel('x')
plt.ylabel('y')
plt.savefig('Fig 1.6.png')
