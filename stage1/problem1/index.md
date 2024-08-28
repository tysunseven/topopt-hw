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

$$c=\sum_{i=1}^{\rm nely}\sum_{j=1}^{\rm nelx}x(i,j)^{\rm penal}U_e(i,j)^TK_eU_e(i,j).$$

- **Constraint**: The problem imposes a volume constraint, ensuring that the material used in the optimized structure does not exceed the initially specified volume fraction. This constraint is expressed as:

$$\sum_{i=1}^{\rm nely}\sum_{j=1}^{\rm nelx}x(i,j)<{\rm volfrac}\times {\rm nelx}\times {\rm nely}.$$




 
- 
This is the report for Problem 1 of Stage 1. In this report, I will analyze the influence of different parameters on the topology of the MBB-beam using the provided Matlab code. The parameters I will investigate include the filter type, filter size, penalization power, and discretization. The results are discussed in the following sections.

## Filter Type Influence
Filter type (ft) can significantly impact the resulting topology. In this study, we have two types of filters: sensitivity filter (ft = 1) and density filter (ft = 2).

<p align="center">
  <img alt="Light" src="../../image/stage1/problem1/MBB_top_nelx60_nely20_rmin1.5_penal3.0_ft1_loop57.png" width="45%">
&nbsp; &nbsp; &nbsp; &nbsp;
  <img alt="Dark" src="../../image/stage1/problem1/MBB_top_nelx60_nely20_rmin1.5_penal3.0_ft2_loop90.png" width="45%">
</p>

In comparing these images, we can observe that the sensitivity filter tends to produce more distinct structural features, whereas the density filter often results in smoother transitions between solid and void regions. This difference occurs because the sensitivity filter directly influences the gradient calculation, leading to sharper boundaries, while the density filter averages densities, resulting in more gradual changes.

Analysis: The choice between sensitivity and density filters depends on the desired balance between sharpness and smoothness in the topology. For designs requiring well-defined features, the sensitivity filter might be preferred. For smoother distributions, the density filter is more suitable.

## Filter Size Influence
The filter size (rmin) determines the extent of the influence of neighboring elements in the optimization process. Larger filter sizes generally lead to smoother and more connected designs, while smaller filter sizes may allow for finer details.

Image: MBB_top_nelx60_nely20_rmin1.5_penal3.0_ft2_loopXX.png vs. MBB_top_nelx60_nely20_rmin3.5_penal3.0_ft2_loopXX.png
When comparing these images, it is evident that a smaller filter size of 1.5 allows for more intricate features, potentially leading to a more optimized but complex structure. Conversely, a larger filter size of 3.5 produces a simpler, more robust structure with fewer fine details.

Analysis: The filter size must be chosen based on the specific application. For complex, lightweight designs, a smaller rmin might be advantageous. For more robust and manufacturable designs, a larger rmin is preferable.

## Penalization Power Influence
The penalization power (penal) influences the material distribution within the design. Higher penal values drive the design towards a more binary (black and white) distribution, which is essential for practical manufacturing.

Image: MBB_top_nelx60_nely20_rmin2.5_penal2.0_ft1_loopXX.png vs. MBB_top_nelx60_nely20_rmin2.5_penal4.0_ft1_loopXX.png
Comparing these images shows that higher penalization (penal = 4.0) results in a more distinct separation between material and void, whereas a lower penalization (penal = 2.0) leads to intermediate grey areas in the topology.

Analysis: For manufacturing purposes, a higher penalization is generally preferred as it ensures a clear material distribution. However, this can lead to less optimal designs if the penalization is too high.

## Discretization Influence
Discretization (nelx * nely) defines the resolution of the finite element mesh. Higher discretization allows for finer details but at the cost of increased computational time.

Image: MBB_top_nelx60_nely20_rmin2.5_penal3.0_ft1_loopXX.png vs. MBB_top_nelx120_nely40_rmin2.5_penal3.0_ft1_loopXX.png
The image with higher discretization (nelx = 120) reveals more detailed features compared to the lower discretization (nelx = 60). The increased resolution provides a more refined structure but requires significantly more computational resources.

Analysis: The choice of discretization depends on the balance between computational efficiency and the level of detail required in the final design. For high-fidelity designs, higher discretization is necessary.

## Analysis of Oscillation Phenomenon
During the study, it was observed that when rmin = 2.5 and 3.5 with penal = 4.0 and nelx = 60, the optimization process failed to converge and resulted in oscillations.

Image: MBB_top_nelx60_nely20_rmin2.5_penal4.0_ft1_error.png and MBB_top_nelx60_nely20_rmin3.5_penal4.0_ft1_error.png
These oscillations are likely caused by a combination of the high penalization power and the selected filter sizes. The high penal value forces the design towards a binary solution, but the filter size and resolution may not adequately support such a binary transition, leading to instability in the optimization process.

Analysis: To avoid such oscillations, either the penalization power should be reduced, or the filter size should be adjusted. Additionally, increasing the discretization might help stabilize the optimization process.

Conclusion
Through this analysis, the influences of different parameters on the topology of the MBB-beam have been explored. Each parameter plays a crucial role in determining the final design, and careful consideration is required to achieve the desired outcome.


