import cmath 
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

class Complex:
    def __init__(self, real, imag=0.0):
        self.real = real
        self.imag = imag

    def __repr__(self):
        if self.imag >= 0:
            return f"{self.real} + {self.imag}i"
        else:
            return f"{self.real} - {abs(self.imag)}i"

def f(z, c):
    return z**3 - 2*z + c

def f_defin(z):
    return 3*z**2 - 2

def newtons_method(c, initial_guess=0.0, tolerance=0.000001, max_iterations=1000):
    z = initial_guess
    for _ in range(max_iterations):
        f_z = f(z, c)
        f_defin_z = f_defin(z)
        z_next = z - f_z / f_defin_z
        
        if abs(z_next - z) < tolerance:  
            return z_next
        z = z_next
    return None

def is_near_duplicate(zero, zeros, tolerance=1e-6):
    for existing_zero in zeros:
        if abs(zero - existing_zero) < tolerance:
            return True
    return False

def julia_set(c, xlim=(-2, 2), ylim=(-2, 2), width=400, height=400, max_iter=50):
    x = np.linspace(xlim[0], xlim[1], width)
    y = np.linspace(ylim[0], ylim[1], height)
    
    Z = x[:, None] + 1j * y[None, :]
    img = np.zeros(Z.shape, dtype=int)

    mask = np.ones(Z.shape, dtype=bool)
    
    for i in range(max_iter):
        Z[mask] = Z[mask] ** 2 + c
        mask = mask & (np.abs(Z) < 2)
        img[mask] = i

    return img

def update_plot(val):
    real_part = real_slider.val
    imag_part = imag_slider.val
    c = complex(real_part, imag_part)

    if val == 0 or c != update_plot.last_c:
        update_plot.last_c = c
        julia_img = julia_set(c)
        julia_ax.clear()
        julia_ax.imshow(julia_img, extent=(-2, 2, -2, 2), cmap='twilight', origin='lower')
        julia_ax.set_title(f'Julia Set for c = {c}')

        real_zeros.clear()
        complex_zeros.clear()
        for guess in initial_guesses:
            zero = newtons_method(c, initial_guess=guess, tolerance=tolerance)
            if zero is not None:
                if abs(zero.imag) < 1e-6:  
                    real_zero = round(zero.real, 6)
                    real_zeros.add(real_zero)
                else:
                    if not is_near_duplicate(zero, complex_zeros):
                        complex_zeros.add(zero)

        newton_ax.clear()
        all_zeros = list(real_zeros) + list(complex_zeros)
        newton_ax.scatter([z.real for z in all_zeros], [z.imag for z in all_zeros], c='red')
        
        newton_ax.set(xlim=(-5, 5), ylim=(-5, 5), aspect='equal')
        newton_ax.spines['bottom'].set_position('zero')
        newton_ax.spines['left'].set_position('zero')
        newton_ax.spines['top'].set_visible(False)
        newton_ax.spines['right'].set_visible(False)
        newton_ax.set_xlabel('x', size=14, labelpad=-24, x=1.03)
        newton_ax.set_ylabel('y', size=14, labelpad=-21, y=1.02, rotation=0)
        newton_ax.grid(which='major', color='grey', linewidth=1, linestyle='-', alpha=0.2)
        newton_ax.legend(['Zeros'])

        plt.draw()
        
update_plot.last_c = complex(0, 0)

def main():
    global real_slider, imag_slider, julia_ax, newton_ax, real_zeros, complex_zeros, initial_guesses, tolerance
    
    tolerance = 0.000001
    initial_guesses = [0.0, 1.0, -1.0, 2.0, -2.0, 1+1j, -1-1j, 2-1j, -2+1j]  
    real_zeros = set()  
    complex_zeros = set()

    fig, (julia_ax, newton_ax) = plt.subplots(1, 2, figsize=(15, 10))

    axcolor = 'lightgoldenrodyellow'
    ax_real = plt.axes([0.2, 0.01, 0.65, 0.03], facecolor=axcolor)
    ax_imag = plt.axes([0.2, 0.05, 0.65, 0.03], facecolor=axcolor)

    real_slider = Slider(ax_real, 'Real (c)', -2.0, 2.0, valinit=0.0)
    imag_slider = Slider(ax_imag, 'Imaginary (c)', -2.0, 2.0, valinit=0.0)

    update_plot(0)

    real_slider.on_changed(update_plot)
    imag_slider.on_changed(update_plot)

    plt.title("f(z) = z^3 - 2*z + c")
    plt.show()

if __name__ == "__main__":
    main()
