# Proof and Comparison: 4f System vs. Single-Lens System

This note derives the field propagation of a classical 4f optical system and a single-lens spatial-filtering system, and clarifies when the 4f configuration implements an ideal Fourier transform while the single-lens system exhibits an unavoidable residual quadratic phase.

---

## 1. 4f System

### 1.1 Propagation from Object Plane to former Lens 1


$$
u_1(x_1,y_1)=
\frac{1}{i\lambda}\frac{e^{ikf}}{f}
\exp\!\left(\frac{ik(x_1^2+y_1^2)}{2f}\right)
\iint_{-\frac{D_0}{2}}^{\frac{D_0}{2}}
u_0(x_0,y_0)
\exp\!\left(\frac{ik(x_0^2+y_0^2)}{2f}\right)
\exp\!\left(\frac{ik(x_0x_1+y_0y_1)}{f}\right)
dx_0\,dy_0.
$$


Multiply by the lens phase

$$
\exp\!\left(-\frac{ik(x_1^2+y_1^2)}{2f}\right).
$$


### 1.2 Propagation from Lens 1 to the Fourier Plane


$$
u_2(x_2,y_2)\propto
\exp\!\left(\frac{ik(x_2^2+y_2^2)}{2f}\right)
\iint_{-\frac{D_1}{2}}^{\frac{D_1}{2}}
\exp\!\left(\frac{ik(x_1^2+y_1^2)}{2f}\right)
\exp\!\left(-\frac{ik(x_1x_2+y_1y_2)}{f}\right)
\left[
\iint_{-\frac{D_0}{2}}^{\frac{D_0}{2}}
u_0(x_0,y_0)
\exp\!\left(\frac{ik(x_0^2+y_0^2)}{2f}\right)
\exp\!\left(-\frac{ik(x_0x_1+y_0y_1)}{f}\right)
dx_0\,dy_0
\right]
dx_1\,dy_1.
$$


Let

$$
g(x_0,y_0)=u_0(x_0,y_0)
\exp\!\left(\frac{ik(x_0^2+y_0^2)}{2f}\right).
$$


---

### 1.3 Method I (Direct Convolution View)


$$
u_2(x_2,y_2)\propto
\exp\!\left(\frac{ik(x_2^2+y_2^2)}{2f}\right)
\iint_{-\frac{D_1}{2}}^{\frac{D_1}{2}}
\exp\!\left(\frac{ik(x_1^2+y_1^2)}{2f}\right)
\exp\!\left(-\frac{ik(x_1x_2+y_1y_2)}{f}\right)
G\!\left(\frac{kx_1}{f},\frac{ky_1}{f}\right)
dx_1\,dy_1.
$$


Formally,

$$
u_2(x_2,y_2)=
\exp\!\left(\frac{ik(x_2^2+y_2^2)}{2f}\right)
\big[g(-x_2,-y_2)\ast
\mathcal F\{\exp(\frac{ik(x_1^2+y_1^2)}{2f})\}\big],
$$

but the kernel \(\mathcal F\{\cdot\}\) is not analytically tractable in this form.

An alternative attempt:

$$
u_2(x_2,y_2)\propto
\iint_{-\frac{D_1}{2}}^{\frac{D_1}{2}}
\exp\!\left(-ik\frac{(x_1-x_2)^2+(y_1-y_2)^2}{f}\right)
G\!\left(\frac{kx_1}{2f},\frac{ky_1}{f}\right)
dx_1\,dy_1
\approx g(-x_2,-y_2)\ast C,
$$

where

$$
C=\iint_{-\frac{D_1}{2}}^{\frac{D_1}{2}}
\exp\!\left(-ik\frac{(x_1-x_2)^2+(y_1-y_2)^2}{f}\right)
dx_1\,dy_1.
$$


If \(D_1\gg D_2\), the bounds extend to infinity. Using

$$
\int_{-\infty}^{\infty} e^{ix^2}\,dx = \sqrt{\frac{\pi}{2}}(1+i),
$$

the kernel becomes a constant, hence this route does not yield a useful spatial filter.

---

### 1.4 Method II (Exchange of Integrals)

Starting again from


$$
u_2(x_2,y_2)\propto
\exp\!\left(\frac{ik(x_2^2+y_2^2)}{2f}\right)
\iint_{-\frac{D_1}{2}}^{\frac{D_1}{2}}
\exp\!\left(\frac{ik(x_1^2+y_1^2)}{2f}\right)
\exp\!\left(-\frac{ik(x_1x_2+y_1y_2)}{f}\right)
\left[\cdots\right]
dx_1\,dy_1,
$$

exchange the order of integration (integrate over \(x_1,y_1\) first):


$$
u_2(x_2,y_2)\propto
\iint_{-\frac{D_0}{2}}^{\frac{D_0}{2}}
u_0(x_0,y_0)
\exp\!\left(\frac{ik(x_0^2+y_0^2)}{2f}\right)
\left[
\iint_{-\frac{D_1}{2}}^{\frac{D_1}{2}}
\exp\!\left(\frac{ik(x_1^2+x_2^2-2x_0x_1-2x_1x_2)}{2f}\right)
\exp\!\left(\frac{ik(y_1^2+y_2^2-2y_0y_1-2y_1y_2)}{2f}\right)
dx_1\,dy_1
\right]
dx_0\,dy_0.
$$


Complete the square to obtain \((x_1-x_0-x_2)^2\) and \((y_1-y_0-y_2)^2\), introducing a term \(-2x_0x_2\) (and similarly for \(y\)):


$$
u_2(x_2,y_2)\propto
\iint_{-\frac{D_0}{2}}^{\frac{D_0}{2}}
u_0(x_0,y_0)
\exp\!\left(-\frac{ik(x_0x_2+y_0y_2)}{f}\right)
\left[
\iint_{-\frac{D_1}{2}}^{\frac{D_1}{2}}
\exp\!\left(\frac{ik(x_1-x_0-x_2)^2}{2f}\right)
\exp\!\left(\frac{ik(y_1-y_0-y_2)^2}{2f}\right)
dx_1\,dy_1
\right]
dx_0\,dy_0.
$$


When \(D_1\gg D_0+D_2\), the inner integral is effectively constant (independent of \(x_0+x_2\)), again using

$$
\int_{-\infty}^{\infty} e^{ix^2}\,dx = \sqrt{\frac{\pi}{2}}(1+i).
$$


Thus

$$
u_2(x_2,y_2)\propto
\iint_{-\frac{D_0}{2}}^{\frac{D_0}{2}}
u_0(x_0,y_0)
\exp\!\left(-\frac{ik(x_0x_2+y_0y_2)}{f}\right)
dx_0\,dy_0
=
U_0\!\left(\frac{kx_2}{f},\frac{ky_2}{f}\right).
$$


**Conclusion:** a 4f system implements an ideal Fourier transform.

---

## 2. Single-Lens System

### 2.1 Propagation from Object Plane to Lens


$$
u_1(x_1,y_1)=
\frac{1}{i\lambda}\frac{e^{iku}}{u}
\exp\!\left(\frac{ik(x_1^2+y_1^2)}{2u}\right)
\iint_{-\frac{D_0}{2}}^{\frac{D_0}{2}}
u_0(x_0,y_0)
\exp\!\left(\frac{ik(x_0^2+y_0^2)}{2u}\right)
\exp\!\left(\frac{ik(x_0x_1+y_0y_1)}{u}\right)
dx_0\,dy_0.
$$


Multiply by the lens phase

$$
\exp\!\left(-\frac{ik(x_1^2+y_1^2)}{2f}\right).
$$


### 2.2 Propagation from Lens to Image Plane


$$
u_2(x_2,y_2)\propto
\exp\!\left(\frac{ik(x_2^2+y_2^2)}{2v}\right)
\iint_{-\frac{D_1}{2}}^{\frac{D_1}{2}}
\exp\!\left[
\frac{ik(x_1^2+y_1^2)}{2}
\left(\frac{1}{u}+\frac{1}{v}-\frac{1}{f}\right)
\right]
\exp\!\left(-\frac{ik(x_1x_2+y_1y_2)}{v}\right)
G\!\left(\frac{kx_1}{u},\frac{ky_1}{u}\right)
dx_1\,dy_1.
$$


For imaging, where u is object distance and v is image distance

$$
\frac{1}{u}+\frac{1}{v}=\frac{1}{f},
$$

and define

$$
g(x_0,y_0)=u_0(x_0,y_0)
\exp\!\left(\frac{ik(x_0^2+y_0^2)}{2u}\right).
$$


Then


$$
u_2(x_2,y_2)\propto
\exp\!\left(\frac{ik(x_2^2+y_2^2)}{2v}\right)
\iint_{-\frac{D_1}{2}}^{\frac{D_1}{2}}
\exp\!\left(-\frac{ik(x_1x_2+y_1y_2)}{v}\right)
G\!\left(\frac{kx_1}{u},\frac{ky_1}{u}\right)
dx_1\,dy_1.
$$


Rewrite the exponential in spatial-frequency variables:


$$
u_2(x_2,y_2)\propto
\exp\!\left(\frac{ik(x_2^2+y_2^2)}{2v}\right)
\iint_{-\frac{D_1}{2}}^{\frac{D_1}{2}}
\exp\!\left[
-i\left(
\frac{kx_1}{u}\frac{ux_2}{v}+
\frac{ky_1}{u}\frac{uy_2}{v}
\right)
\right]
G\!\left(\frac{kx_1}{u},\frac{ky_1}{u}\right)
dx_1\,dy_1
=
\exp\!\left(\frac{ik(x_2^2+y_2^2)}{2v}\right)
g\!\left(-\frac{ux_2}{v},-\frac{uy_2}{v}\right).
$$


**Conclusion:** the single-lens system contains an inevitable residual quadratic phase

$$
\exp\!\left(\frac{ik(x_2^2+y_2^2)}{2v}\right),
$$

which represents a fundamental imaging imperfection compared with the 4f case.

---
