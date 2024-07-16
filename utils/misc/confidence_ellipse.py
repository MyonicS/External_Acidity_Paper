from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
import numpy as np



def confidence_ellipse_matrix(x_mean,y_mean,matrix,ax, facecolor='none',n_std=1, **kwargs):
    cov = matrix
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensional dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                      facecolor=facecolor, **kwargs)

    # Calculating the standard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.

    scale_x = np.sqrt(cov[0, 0]) * n_std/1000*8.314

    # calculating the standard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(x_mean, y_mean)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)
