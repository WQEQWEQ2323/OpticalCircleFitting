The main script is **`four_f_F100mm_formula.m`**; all other `.m` files are utility/tool functions.

1. The word *formula* in the filename indicates that this script uses the **first-order Hankel transform** of the sombrero function.(transmission_function.m)
   Similar to the fact that the **zero-order Hankel transform** of a sombrero function yields an Airy function, the **first-order Hankel transform** of the sombrero function is
   $$
   \mathrm{result}
   =\frac{\pi R}{2}\,
   \frac{
     J_1(rR)\,H_0(rR)-J_0(rR)\,H_1(rR)
   }{r},
   $$
   
   where \(r\) and \(\phi\) are the polar coordinates on the image plane:
   $$
   \rho = \sqrt {x^2 + y^2},\qquad
   \phi = \arctan\!\left(\frac{y}{x}\right),
   $$
   and \(k\) is the wave vector. Here \(J_0(\cdot)\) and \(J_1(\cdot)\) denote the zero- and first-order Bessel functions, and \(H_0(\cdot)\) and \(H_1(\cdot)\) denote the zero- and first-order Struve functions, respectively. The wavenumber is
   $$
   k=\frac{2\pi}{\lambda}.
   $$

**Functions of the other `.m` files:**

- **`F_ZerosFilling.m`** and **`F_UnZerosFilling.m`**  
  These functions are used for **zero-padding** and **matrix cropping**, respectively.

- **`jfft2`** and **`jifft2`**  
  These are used to perform **angular spectrum method (ASM)** simulations for optical field propagation.

- **Struve/Bessel-related functions**  
  Functions such as those computing **Struve functions** are utility routines for constructing and evaluating the **transmission function**.