<script type="text/javascript" async
  src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js">
</script>
# Problem 1: Test influence of discretization, filter type, filter size and penalization power
- **Design domain**: A rectangle with an aspect ratio of 3:1, representing half of an MBB beam.
- **Boundary conditions**: The x-component on the left edge is fixed, and the y-component at the bottom right corner is fixed.
```matlab
fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
```
- **Loading conditions**: A vertical downward force is applied at the top-left corner.
```matlab
F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1);
```
## Optimization Model: SIMP
The SIMP (Solid Isotropic Material with Penalization) method is a widely used approach in topology optimization. It defines the material properties as a function of the design variables, typically representing material density. The method penalizes intermediate densities to push the design towards a clear distinction between solid and void regions, enabling efficient and manufacturable structures. The primary goal of SIMP is to optimize material distribution within a given design domain to achieve maximum stiffness under specified constraints.
- **Design Variable**: Material density `x`, where `x` ranges from 0 to 1, with 0 indicating void and 1 indicating solid material. In the code, `x` is implemented as an `nely` by `nelx` matrix, where each element represents the density of a specific element in the design domain. The values of this matrix are initially set to `volfrac`, the prescribed volume fraction, to ensure a uniform distribution of material density at the start of the optimization process.
- **Objective Function**: Minimize compliance (maximize stiffness) of the structure, represented as

$$c(x)=U^TKU=\sum_{e=1}^{\rm nelx\times nely}E_e(x_e)U_e^TK_0U_e,\quad E_e(x_e)=\begin{cases}x_e^p&\text{in 99-line code}\\E_{\rm min}+x_e^p(E_0-E_{\rm min})&\text{in 88-line code}\end{cases}$$
- **Constraint**: The problem imposes a volume constraint, ensuring that the material used in the optimized structure does not exceed the initially specified volume fraction. This constraint is expressed as:

$$\sum_{e=1}^{\rm nely\times nelx}x_e<{\rm volfrac}\times {\rm nelx}\times {\rm nely}.$$

## Sensitivity Analysis
In the context of topology optimization, sensitivity analysis is essentially the process of calculating the derivative of the objective function with respect to the design variables. This derivative, or sensitivity, provides the necessary information to guide the optimization process, ensuring that updates to the design variables are made in a direction that drives the solution toward an optimal configuration.

$$\begin{align*}
\frac{\partial c}{\partial x_e}&=\frac{\partial U^T}{\partial x_e}KU+U^T\frac{\partial K}{\partial x_e}U+U^TK\frac{\partial U}{\partial x_e}=U^T\frac{\partial K}{\partial x_e}U+2U^TK\frac{\partial U}{\partial x_e} \\
&=U^T\frac{\partial K}{\partial x_e}U-2U^TKK^{-1}\frac{\partial K}{\partial x_e}U=-U^T\frac{\partial K}{\partial x_e}U=-p x_e^{p-1}（E_0-E_{\min}） U_e^T K_0 U_e
\end{align*}$$

When I first saw \\(-p x_e^{p-1}(E_0-E_{\min})\\), my initial reaction was that this term comes from the derivative of \\(E_e(x_e)\\) with respect to \\(x_e\\), with an extra negative sign added for some reason unknown to me. As a result, I mistakenly interpreted the \\(\text{dc}\\) in the code as the negative of the derivative of \\(c\\) with respect to \\(x\\). However, in reality, \\(U_e\\) also contains the dependency of \\(c\\) on \\(x\\), and when all these dependencies are taken into account, the correct result just happens to introduce a negative sign. This means that \\(\text{dc}\\) is indeed the derivative of \\(c\\) with respect to \\(x\\), not the negative of the derivative as I initially understood.

We can see that derivatives of compliance with respect to every position are all negative. This property reflects the fact that increase of element densities always leads to less compliant structure.



## Optimization Algorithm: Optimality Criteria (OC)

The Optimality Criteria (OC) method is the core optimization algorithm used in this code to iteratively adjust the design variables, guiding the solution towards optimality. In each iteration, the OC method updates the design variables x to minimize the objective function (compliance) while ensuring the volume constraint is met. This process is governed by the Lagrangian function, which integrates the objective with the optimization constraints.

 $$\mathcal{L} = c + \lambda\left(\sum_{e}^{\rm nely\times nelx} x_e - {\rm volfrac}\times {\rm nelx}\times {\rm nely}\right) + \sum_{e} \lambda_e^+ (x_e - 1) + \sum_{e} \lambda_e^- (x_{\text{min}} - x_e).$$

The Karush-Kuhn-Tucker (KKT) Condition of this Lagrangian function reads

$$ \frac{\partial c}{\partial x_e}+\lambda+\sum_e\lambda_e^+-\sum_e\lambda_e^-=0 $$

The idea of the Optimality Criteria (OC) method is to approximate \\(c\\) for a given initial guess solution \\(x_0\\) by

$$ \tilde{c}_{x_0}(x)=c(x_0)-\sum_{e}\frac{\partial c}{\partial x_e}(x_0)x^2_{e0}\left(\frac{1}{x_e}-\frac{1}{x_{e0}}\right),\quad \frac{\partial \tilde{c}_{x_0}}{\partial x_e}=\frac{\partial c}{\partial x_e}(x_0)\frac{x^2_{e_0}}{x_e^2} $$

Then we have new KKT condition for new Lagrangian reads

$$ \frac{\partial c}{\partial x_e}(x_0)\frac{x^2_{e_0}}{x_e^2}+\lambda+\sum_e\lambda_e^+-\sum_e\lambda_e^-=0\Longrightarrow x_e=\sqrt{\frac{B_e}{\tilde{\lambda}}}x_{e0},\quad B_e=-\frac{\partial c}{\partial x_e}(x_0),\quad \tilde{\lambda}=\lambda+\sum_e\lambda_e^+-\sum_e\lambda_e^- $$

This way, we obtain a better solution \\(x_1\\) from the initial guess \\(x_0\\), and then \\(x_1\\) in turn becomes the new guess for the next iteration, thus starting the iterative process. In practice, \\(x_0\\) is typically chosen as a uniform distribution, and we constrain each iteration step so that the change does not exceed a specified value \\( \text{move} \\), while also ensuring that the result remains within the bounds of 0 and 1. This leads to the code implementation of

$$ x_{\text{new}} = \max\left(0, \max\left(x - \text{move}, \min\left(1, \min\left(x + \text{move}, x \cdot \sqrt{\frac{-dc}{dv \cdot \lambda}}\right)\right)\right)\right) $$

In practice, we also need to apply an additional filtering step to \\(dc\\) and \\(dv\\) to avoid numerical instability and mesh dependency issues, which will be discussed in detail in later sections. We use the bisection method to find the value of \\(\lambda\\) to satisfy the volume constraint. The iterative process continues until the maximum change in the design variables between iterations falls below a specified threshold (e.g., 0.01). This indicates that the solution has converged to a near-optimal design.

We haven't provided a formal proof that \\(x_1\\) is better than \\(x_0\\), but I can offer some intuitive explanations. The total volume remains constant before and after the iteration, so the density increases in some positions and decreases in others. The total increase in density equals the total decrease. The factor \\(\lambda\\) serves as a common coefficient and can be disregarded for now. The degree of increase or decrease in each position is directly proportional to the square root of \\(-\text{dc}\\), which in turn is proportional to the absolute value of the derivative of the objective function with respect to that position. Thus, positions where the material density increases have a higher absolute value of the derivative compared to positions where the density decreases. According to sensitivity analysis, this means the objective function is more sensitive to changes in these positions. Therefore, even though the total increase in density equals the total decrease, the effect on the objective function from the increase is greater than that from the decrease. Since an increase in material density corresponds to a decrease in the objective function, the value of the objective function will be lower after the iteration.

## Filter type: Sensitivity Filtering
In topology optimization, filtering plays a critical role in ensuring the stability and physical meaningfulness of the optimization process. Without filtering, the optimization algorithm can produce results that exhibit numerical instabilities, such as the formation of "checkerboard" patterns or other non-physical artifacts. These issues arise because the optimization can exploit the high sensitivity of certain elements to create designs that are mathematically optimal but physically unrealistic. Filtering smooths out these variations by considering the influence of neighboring elements, thereby promoting a more stable and reliable convergence toward a manufacturable design.

In Sensitivity Filtering, we first define a sparse matrix \\( H \\), which provides the weights for performing a weighted average of the values at a point with the values of points in a nearby region selected in a certain manner. The method for selecting the region and the weighting of the values can be adjusted according to your specific needs.

In Sensitivity Filtering, the object being weighted is the sensitivity \\( dc \\). First, \\( dc \\) is multiplied by \\( x \\) at each position, and the resulting product is then multiplied by the matrix \\( H \\). This result is then divided by the normalization factor \\( Hs \\), and finally, each position is divided by \\( x \\). Thus, we effectively perform two layers of weighting: the first layer reduces the influence of positions with low material density, and the second layer averages each point with the surrounding positions.

In practice, we divide each position by the larger value between \\(1e-3\\) and \\(x\\) to avoid division by zero.
```matlab
if ft == 1
  dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
end
```



## Filter type: Density Filtering
Compared to Sensitivity Filtering, Density Filtering simplifies the processing of \\(dc\\) by removing one layer of weighting by \\(x\\), using only the matrix \\(H\\) to weight \\(dc\\). However, Density Filtering introduces an additional step: when using the bisection method to find the value of \\(\lambda\\), it compares the volume of \\(x_{\text{Phys}}\\) (the weighted result of \\(x_{\text{new}}\\) by \\(H\\)) with the volume constraint, rather than comparing the volume of \\(x_{\text{new}}\\) calculated using the previous \\(\lambda\\). Nonetheless, the value used to update \\(x\\) is still \\(x_{\text{new}}\\), not the weighted average \\(x_{\text{new}}\\). The purpose of this operation remains somewhat mysterious to me.
