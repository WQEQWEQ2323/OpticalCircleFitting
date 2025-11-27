# Fresnel Diffraction, Hankel Transform, and Spiral Phase Contrast Imaging

This note derives the Fresnel diffraction imaging process, reviews the 2D Fourier transform in Cartesian coordinates, introduces the n-th order Hankel transform, and shows how it leads to a polar-form Fourier representation. Finally, a compact expression for spiral phase contrast imaging is given.

---

## 1. Fresnel Diffraction Imaging

The Fresnel diffraction integral from the object plane to the observation plane is


$$
u_1(x_1,y_1)=-\frac{i}{\lambda}\iint_{-\infty}^{\infty}u_0(x_0,y_0)\frac{e^{ikR}}{R}\,dx_0\,dy_0,
$$

where

$$
R^2=(x_1-x_0)^2+(y_1-y_0)^2+u^2=\rho^2+u^2,
\quad
\rho^2=(x_1-x_0)^2+(y_1-y_0)^2.
$$


Using the Fresnel approximation,

$$
R=\sqrt{\rho^2+u^2}
= u\sqrt{\left(\frac{\rho}{u}\right)^2+1}
\approx u+\frac{\rho^2}{2u}.
$$


Define the field just after multiplying by the lens quadratic phase as

$$
\hat u_1(x_1,y_1)=u_1(x_1,y_1)\exp\!\left(-\frac{ik(x_1^2+y_1^2)}{2f}\right).
$$


Then

$$
u_1(x_1,y_1)=-\frac{ie^{iku}}{\lambda u}
\iint_{-\infty}^{\infty}u_0(x_0,y_0)
\exp\!\left(\frac{ik\rho^2}{2u}\right)\,dx_0\,dy_0.
$$


Expanding \(\rho^2\):

$$
\rho^2=x_1^2+y_1^2+x_0^2+y_0^2-2x_1x_0-2y_1y_0,
$$

gives

$$
u_1(x_1,y_1)=
-\frac{ie^{iku}}{\lambda u}
\exp\!\left(\frac{ik(x_1^2+y_1^2)}{2u}\right)
\iint_{-\infty}^{\infty}
u_0(x_0,y_0)
\exp\!\left(\frac{ik(x_0^2+y_0^2)}{2u}\right)
\exp\!\left[-i2\pi\!\left(x_0\frac{x_1}{\lambda u}+y_0\frac{y_1}{\lambda u}\right)\right]
dx_0\,dy_0.
$$


Let

$$
g(x_0,y_0)=u_0(x_0,y_0)
\exp\!\left(\frac{ik(x_0^2+y_0^2)}{2u}\right),
$$

and denote its Fourier transform as \(G(\cdot,\cdot)\). Then

$$
u_1(x_1,y_1)=
-\frac{ie^{iku}}{\lambda u}
\exp\!\left(\frac{ik(x_1^2+y_1^2)}{2u}\right)
G\!\left(\frac{x_1}{\lambda u},\frac{y_1}{\lambda u}\right).
$$


After another propagation distance \(v\), the field becomes


$$
u_2(x_2,y_2)=
-\frac{ie^{ikv}}{\lambda v}
\exp\!\left(\frac{ik(x_2^2+y_2^2)}{2v}\right)
\iint_{-\infty}^{\infty}
\hat u_1(x_1,y_1)
\exp\!\left(\frac{ik(x_1^2+y_1^2)}{2v}\right)
G\!\left(\frac{x_1}{\lambda u},\frac{y_1}{\lambda u}\right)
\exp\!\left[-i2\pi\!\left(x_1\frac{x_2}{\lambda v}+y_1\frac{y_2}{\lambda v}\right)\right]
dx_1\,dy_1.
$$


With the imaging condition

$$
\frac{1}{u}+\frac{1}{v}=\frac{1}{f},
$$

one obtains the standard scaled imaging relation


$$
u_2(x_2,y_2)=
-\frac{u}{v}e^{ik(u+v)}
\exp\!\left(\frac{ik(x_2^2+y_2^2)}{2v}\right)
g\!\left(-\frac{u}{v}x_2,-\frac{u}{v}y_2\right).
$$


---

## 2. 2D Fourier Transform in Cartesian Coordinates


$$
F(\xi,\eta)=\iint_{-\infty}^{\infty}
f(x,y)e^{-i2\pi(x\xi+y\eta)}\,dx\,dy,
$$


$$
f(x,y)=\iint_{-\infty}^{\infty}
F(\xi,\eta)e^{i2\pi(\xi x+\eta y)}\,d\xi\,d\eta.
$$


Fourier pair notation:

$$
f(x,y)\leftrightarrow F(\xi,\eta),
\qquad
F(\xi,\eta)\leftrightarrow f(-x,-y).
$$


---

## 3. Hankel Transform

Using the Bessel function identity

$$
J_n(x)=\frac{1}{2\pi}\int_{-\pi}^{\pi}e^{i(n\theta-x\sin\theta)}\,d\theta,
$$

the n-th order Hankel transform is defined as


$$
F(\rho)=\mathcal H_n\{f(r)\}
=\int_{0}^{\infty} f(r)\,r\,J_n(\rho r)\,dr,
$$

with inverse

$$
f(r)=2\pi\int_{0}^{\infty}
F(\rho)\,\rho\,J_n(2\pi\rho r)\,d\rho.
$$


---

## 4. 2D Fourier Transform in Polar Coordinates

Start from the Fourier transform

$$
F(\xi,\eta)=\iint_{-\infty}^{\infty} f(x,y)\,e^{-i(x\xi+y\eta)}\,dx\,dy.
$$


Let

$$
x=r\cos\theta,\quad y=r\sin\theta,
\qquad
\xi=\rho\cos\phi,\quad \eta=\rho\sin\phi.
$$


Then

$$
F(\rho,\phi)=\int_{0}^{2\pi}\!d\theta\int_{0}^{\infty}
f(r,\theta)\,r\,e^{-i\rho r\cos(\theta-\phi)}\,dr.
$$


Expand \(f(r,\theta)\) in angular harmonics:

$$
f(r,\theta)=\sum_{n=-\infty}^{\infty} f_n(r)e^{in\theta},
\quad
f_n(r)=\frac{1}{2\pi}\int_{0}^{2\pi}f(r,\theta)e^{-in\theta}\,d\theta.
$$


Substituting yields

$$
F(\rho,\phi)=
2\pi\sum_{n=-\infty}^{\infty}
i^{-n}e^{in\phi}
\int_{0}^{\infty} f_n(r)\,r\,J_n(\rho r)\,dr.
$$


Define

$$
F_n(\rho)=2\pi i^{-n}\mathcal H_n\{f_n(r)\},
$$

so that

$$
F(\rho,\phi)=\sum_{n=-\infty}^{\infty}F_n(\rho)e^{in\phi}.
$$


Conversely,

$$
f_n(r)=\frac{i^{n}}{2\pi}\mathcal H_n^{-1}\{F_n(\rho)\}.
$$


---

## 5. Spiral Phase Contrast Imaging (SPCI)

Consider a circular aperture carrying a unit spiral phase:


$$
f(r,\theta)=\mathrm{circle}\!\left(\frac{r}{R}\right)e^{i\theta},
$$

where

$$
\mathrm{circle}(r)=
\begin{cases}
1, & r>1,\\
0, & 0\le r \le 1.
\end{cases}
$$


Its Fourier transform in polar form is


$$
F(\rho,\phi)=2\pi i^{-1}\,\mathcal H_1\!\left\{\mathrm{circle}\!\left(\frac{r}{R}\right)\right\}.
$$


This compactly summarizes the SPCI transfer behavior in the Hankel domain.
