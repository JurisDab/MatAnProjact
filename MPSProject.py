import cmath 
import numpy as np
import matplotlib.pyplot as plt

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

def main():
    user_input = input("Enter c: ")
    c = float(user_input)
    
    tolerance = 0.000001
    
    initial_guesses = [0.0, 1.0, -1.0, 2.0, -2.0]  
    
    real_zeros = set() 

    for guess in initial_guesses:
        zero = newtons_method(c, initial_guess=guess, tolerance=tolerance) 
        if zero is not None:
            if abs(zero.imag) < 0.000001:
                real_zero = round(zero.real, 6) 
                real_zeros.add(real_zero)  
    
    if real_zeros:
        print("Zeros found:")
        for zero in real_zeros:  
            print(f"Zero: {zero}")
    else:
        print("No zeros found.")

    if real_zeros:
        xmin, xmax, ymin, ymax = -5, 5, -5, 5
        ticks_frequency = 1

        real_zeros_array = np.array(list(real_zeros))
        zeros_array = np.zeros_like(real_zeros_array)
        
        fig, ax = plt.subplots(figsize=(10, 10))
        ax.scatter(real_zeros_array, zeros_array, c='r')

        ax.set(xlim=(xmin-1, xmax+1), ylim=(ymin-1, ymax+1), aspect='equal')

        ax.spines['bottom'].set_position('zero')
        ax.spines['left'].set_position('zero')

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        ax.set_xlabel('x', size=14, labelpad=-24, x=1.03)
        ax.set_ylabel('y', size=14, labelpad=-21, y=1.02, rotation=0)

        x_ticks = np.arange(xmin, xmax+1, ticks_frequency)
        y_ticks = np.arange(ymin, ymax+1, ticks_frequency)
        ax.set_xticks(x_ticks[x_ticks != 0])
        ax.set_yticks(y_ticks[y_ticks != 0])

        ax.grid(which='major', color='grey', linewidth=1, linestyle='-', alpha=0.2)

        plt.show()

if __name__ == "__main__":
    main()
