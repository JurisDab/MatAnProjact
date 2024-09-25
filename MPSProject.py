import cmath 

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

if __name__ == "__main__":
    main()
