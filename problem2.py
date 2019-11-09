import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal

def rgb2bayer(image):
    """Convert image to bayer pattern:
    [B G]
    [G R]

    Args:
        image: Input image as (H,W,3) numpy array

    Returns:
        bayer: numpy array (H,W,3) of the same type as image
        where each color channel retains only the respective 
        values given by the bayer pattern
    """
    assert image.ndim == 3 and image.shape[-1] == 3

    # otherwise, the function is in-place
    bayer = image.copy()

    #
    # You code goes here
    #

    [H,W]=bayer[:,:,0].shape
    bayer_r = np.zeros((H, W))
    bayer_g = np.zeros((H, W))
    bayer_b = np.zeros((H, W))
    for i in range(H):
        for j in range(W):
            if((i%2==1)and(j%2==1)):
                bayer_r[i,j]=bayer[i,j,0]
            elif ((i%2==0)and(j%2==0)):
                bayer_b[i,j]=bayer[i,j,2]
            else:
                bayer_g[i,j]=bayer[i,j,1]

    bayer[:, :, 0] = bayer_r
    bayer[:, :, 1] = bayer_g
    bayer[:, :, 2] = bayer_b

    # plt.imshow(bayer)
    # plt.show()

    assert bayer.ndim == 3 and bayer.shape[-1] == 3
    return bayer

def nearest_up_x2(x):
    """Upsamples a 2D-array by a factor of 2 using nearest-neighbor interpolation.

    Args:
        x: 2D numpy array (H, W)

    Returns:
        y: 2D numpy array if size (2*H, 2*W)
    """
    assert x.ndim == 2
    h, w = x.shape

    #
    # You code goes here
    #
    y = np.empty((2*h, 2*w))

    assert y.ndim == 2 and \
            y.shape[0] == 2*x.shape[0] and \
            y.shape[1] == 2*x.shape[1]
    return y

def bayer2rgb(bayer):
    """Interpolates missing values in the bayer pattern.
    Note, green uses bilinear interpolation; red and blue nearest-neighbour.

    Args:
        bayer: 2D array (H,W,C) of the bayer pattern
    
    Returns:
        image: 2D array (H,W,C) with missing values interpolated
        green_K: 2D array (3, 3) of the interpolation kernel used for green channel
        redblue_K: 2D array (3, 3) using for interpolating red and blue channels
    """
    assert bayer.ndim == 3 and bayer.shape[-1] == 3

    #
    # You code goes here
    #
    image = bayer.copy()
    rb_k = np.array([[0, 1/4, 0], [1/4, 0, 1/4], [0, 1/4, 0]])
    g_k = np.array([[0, 1/4, 0], [1/4, 0, 1/4], [0, 1/4, 0]])
    image_r = signal.convolve2d(image[:, :, 0], rb_k, boundary='symm', mode='same', fillvalue=0)
    image_g = signal.convolve2d(image[:, :, 1], g_k, boundary='symm', mode='same', fillvalue=0)
    image_b = signal.convolve2d(image[:, :, 2], rb_k, boundary='symm', mode='same', fillvalue=0)
    image[:, :, 0] = image_r
    image[:, :, 1] = image_g
    image[:, :, 2] = image_b




    assert image.ndim == 3 and image.shape[-1] == 3 and \
                g_k.shape == (3, 3) and rb_k.shape == (3, 3)
    return image, g_k, rb_k

def scale_and_crop_x2(bayer):
    """Upscamples a 2D bayer pattern by factor 2 and takes the central crop.

    Args:
        bayer: 2D array (H, W) containing bayer pattern

    Returns:
        image_zoom: 2D array (H, W) corresponding to x2 zoomed and interpolated 
        one-channel image
    """
    assert bayer.ndim == 2

    #
    # You code goes here
    #
    cropped = bayer.copy()

    # print(bayer.shape)
    # print(cropped.shape)

    [H, W] = cropped.shape
    scale_up=np.zeros((2*H,2*W))

    for i in range(2*H):
        for j in range(2*W):
            scale_up[i,j]=cropped[i//2,j//2]

    cropped = scale_up[int(H/2):int(3*H/2),int(W/2):int(3*W/2)]

    # plt.imshow(cropped)
    # plt.show()

    assert cropped.ndim == 2
    return cropped
