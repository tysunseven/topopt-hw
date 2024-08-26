# 学习 99 行拓扑优化代码

```matlab
x(1:nely,1:nelx) = volfrac;
```
- 初始化一个矩阵 `x`。这个矩阵表示拓扑优化问题中的设计变量（通常为材料密度），其中每个元素的初始值都设置为 `volfrac`。
- `1:nely` 和 `1:nelx` 是行和列的索引范围，分别表示从第 1 行到第 `nely` 行、从第 1 列到第 `nelx` 列。
- 通过这种索引方式，`x(1:nely, 1:nelx)` 表示选取矩阵 `x` 中从第 1 行到第 `nely` 行、从第 1 列到第 `nelx` 列的所有元素。由于 `x` 是一个尚未初始化的矩阵，这行代码实际上是定义了一个大小为 `nely x nelx` 的矩阵 `x`。
- `= volfrac;`这部分语句将矩阵 `x` 的所有选定元素赋值为 `volfrac`。`volfrac` 是一个标量，表示材料体积分数。通过这行代码，`x` 矩阵的每个元素都被初始化为 `volfrac` 的值。

```matlab
change = 1.;
```
- \1. 表示浮点数 1.0，即使没有写出完整的小数部分（.0），MATLAB 仍然将其视为浮点数。
- change 变量用于控制迭代的停止条件。它初始化为一个非零值（如 1.0）。随后，在每次迭代结束时，change 会被重新计算，以反映当前迭代中设计变量的变化程度。

```matlab
while change > 0.01
  loop = loop + 1;
  xold = x;
  ......
end```

```matlab
% 有限元分析
  [U]=FE(nelx,nely,x,penal);  
```
- 方括号 [] 用于定义函数的输出变量。在这个例子中，U 是 FE 函数的输出变量。因为只有一个输出，U 也是一个单一变量。方括号在这种情况下并不是严格必要的，但在 MATLAB 中，即使只有一个输出变量，使用方括号也是一种常见的习惯。
- 这一行代码的作用是调用 FE 函数，将设计域的尺寸（nelx 和 nely）、设计变量矩阵 x 以及惩罚因子 penal 作为输入，计算结构在当前设计状态下的位移场，并将结果存储在变量 U 中。

```matlab
function [KE]=lk
E = 1.; 
nu = 0.3;
k=[ 1/2-nu/6   1/8+nu/8 -1/4-nu/12 -1/8+3*nu/8 ... 
   -1/4+nu/12 -1/8-nu/8  nu/6       1/8-3*nu/8];
KE = E/(1-nu^2)*[ k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)
                  k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
                  k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
                  k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
                  k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
                  k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
                  k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
                  k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];
```
- 这个函数 lk 计算的是一个二维四节点正方形单元的刚度矩阵（stiffness matrix），该矩阵在有限元分析中用于描述单元的弹性行为。
- 这个函数 lk 没有输入参数，返回一个 8x8 的刚度矩阵 KE，这个矩阵表示单个四节点正方形单元在平面应变条件下的弹性刚度。
- E：定义了材料的弹性模量（Young's modulus）。在这个例子中，E = 1，代表一个单位化的弹性模量。
- nu：定义了材料的泊松比（Poisson's ratio）。泊松比 nu = 0.3 是材料的一种弹性性质，反映了材料在拉伸或压缩时的横向变形与轴向变形的比例。
- k 是一个向量，存储了单元刚度矩阵的某些关键系数，这些系数基于材料属性 E 和 nu 计算得出。
- 这个向量 k 包含了在生成刚度矩阵时所需的8个元素。这些元素是根据材料属性和单元几何形状导出的。
- KE 矩阵的元素由向量 k 中的值填充，排列顺序遵循单元的对称性和力学平衡条件。
- 矩阵 KE 是一个 8x8 的对称矩阵，它的元素描述了单元节点之间的刚度关系。这里的 8 表示四节点单元在二维平面内的自由度总数（每个节点有 2 个自由度：x 和 y）。
- 刚度矩阵 KE 是有限元分析中的核心，表示单个单元的力学响应。这个矩阵乘以节点位移向量，产生单元的内力向量。
- 在整个结构中，所有单元的刚度矩阵会被组合成全局刚度矩阵，最终用于求解结构在外部载荷下的位移场。
```matlab
function [U]=FE(nelx,nely,x,penal)
[KE] = lk; 
K = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1));
F = sparse(2*(nely+1)*(nelx+1),1); U = zeros(2*(nely+1)*(nelx+1),1);
for elx = 1:nelx
  for ely = 1:nely
    n1 = (nely+1)*(elx-1)+ely; 
    n2 = (nely+1)* elx   +ely;
    edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
    K(edof,edof) = K(edof,edof) + x(ely,elx)^penal*KE;
  end
end
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
F(2,1) = -1;
fixeddofs   = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
alldofs     = [1:2*(nely+1)*(nelx+1)];
freedofs    = setdiff(alldofs,fixeddofs);
% SOLVING
U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);      
U(fixeddofs,:)= 0;
```

```matlab
 
% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  [KE] = lk;
  c = 0.;
  for ely = 1:nely
    for elx = 1:nelx
      n1 = (nely+1)*(elx-1)+ely; 
      n2 = (nely+1)* elx   +ely;
      Ue = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],1);
      c = c + x(ely,elx)^penal*Ue'*KE*Ue;
      dc(ely,elx) = -penal*x(ely,elx)^(penal-1)*Ue'*KE*Ue;
    end
  end
% FILTERING OF SENSITIVITIES
  [dc]   = check(nelx,nely,rmin,x,dc);    
% DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
  [x]    = OC(nelx,nely,x,volfrac,dc); 
% PRINT RESULTS
  change = max(max(abs(x-xold)));
  disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',c) ...
       ' Vol.: ' sprintf('%6.3f',sum(sum(x))/(nelx*nely)) ...
        ' ch.: ' sprintf('%6.3f',change )])
% PLOT DENSITIES  
  colormap(gray); imagesc(-x); axis equal; axis tight; axis off;pause(1e-6);
```


