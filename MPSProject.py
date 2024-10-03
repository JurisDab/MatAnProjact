import cmath 
import numpy as np
import matplotlib.pyplot as plt

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

def main():
    user_input = input("Enter c: ")  
    try:
        c = complex(user_input)
    except ValueError:
        print("Invalid input.")
        return

    tolerance = 0.000001
    initial_guesses = [0.0, 1.0, -1.0, 2.0, -2.0, 1+1j, -1-1j, 2-1j, -2+1j]  
    
    real_zeros = set()  
    complex_zeros = set()

    is_c_real = abs(c.imag) < 1e-6

    for guess in initial_guesses:
        zero = newtons_method(c, initial_guess=guess, tolerance=tolerance)
        if zero is not None:
            if abs(zero.imag) < 1e-6:  
                real_zero = round(zero.real, 6)
                real_zeros.add(real_zero)
            else:
                if not is_c_real: 
                    if not is_near_duplicate(zero, complex_zeros):
                        complex_zeros.add(zero)

    if real_zeros:
        print("Zeros found:")
        for zero in real_zeros:  
            print(f"Zero: {zero}")
    
    if complex_zeros:
        print("Zeros found:")
        for zero in complex_zeros:  
            print(f"Zero: {zero}")

    if real_zeros or complex_zeros:
        xmin, xmax, ymin, ymax = -5, 5, -5, 5
        ticks_frequency = 1

        real_zeros_array = np.array(list(real_zeros))
        complex_zeros_real = [z.real for z in complex_zeros]
        complex_zeros_imag = [z.imag for z in complex_zeros]

        fig, ax = plt.subplots(figsize=(10, 10))
        ax.scatter(real_zeros_array, np.zeros_like(real_zeros_array), c='red')
        if complex_zeros:
            ax.scatter(complex_zeros_real, complex_zeros_imag, c='red')

        ax.set(xlim=(xmin-1, xmax+1), ylim=(ymin-1, ymax+1), aspect='equal')
        ax.spines['bottom'].set_position('zero')
        ax.spines['left'].set_position('zero')
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        x_ticks = np.arange(xmin, xmax+1, ticks_frequency)
        y_ticks = np.arange(ymin, ymax+1, ticks_frequency)
        
        ax.set_xlabel('x', size=14, labelpad=-24, x=1.03)
        ax.set_ylabel('y', size=14, labelpad=-21, y=1.02, rotation=0)
        
        ax.set_xticks(x_ticks[x_ticks != 0])
        ax.set_yticks(y_ticks[y_ticks != 0])

        ax.grid(which='major', color='grey', linewidth=1, linestyle='-', alpha=0.2)

        plt.show()

if __name__ == "__main__":
    main()
