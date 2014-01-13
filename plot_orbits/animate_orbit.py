'''
Created on Apr 16, 2013

@author: Tripod
'''


import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation as animation

def animate_orbit(x_pos, y_pos, fig):
    

    line, = plt.plot([], [], 'r-')


# animation function.  This is called sequentially
    def update_line(num, x_pos, y_pos, line):
        line.set_data(x_pos[:num], y_pos[:num])
        return line,

# call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, update_line, fargs=(x_pos, y_pos, line), \
                                   frames=500, interval=1.0, blit=True, repeat=False)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
#anim.save('basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

    return anim