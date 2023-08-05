import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.colors import Normalize, LinearSegmentedColormap

# Ellipse parameters and values
num_ellipses = 10
r = 5  # Radius of the circle containing ellipses
ellipse_params = np.random.rand(num_ellipses, 5)  # h, k, a, b, alpha
values = np.random.rand(num_ellipses)  # Values associated with ellipses

# Create colormap
cmap = plt.cm.get_cmap('viridis')  # Change 'viridis' to any desired colormap

# Normalize values
norm = Normalize(vmin=min(values), vmax=max(values))

# Create figure and axis
fig, ax = plt.subplots()

# Plot circle
circle = plt.Circle((0, 0), r, color=cmap(norm(np.mean(values))))
ax.add_artist(circle)

# Plot ellipses
for params, value in zip(ellipse_params, values):
    h, k, a, b, alpha = params
    ellipse = Ellipse((h, k), 2 * a, 2 * b, angle=np.degrees(alpha), color=cmap(norm(value)))
    ax.add_artist(ellipse)

# Add colorbar
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])  # Set an empty array
cbar = plt.colorbar(sm)
cbar.set_label('Values')

# Set labels and title
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_title('Ellipses within Circle with Colored Shading')

# Set aspect ratio to equal
ax.set_aspect('equal')

# Show plot
plt.show()
