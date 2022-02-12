# Boundary Element Methods By Sauter and Schwab (椭圆边值问题的积分方程)
## 历史背景与引言
* 积分方程方法具有悠久的历史，可以追溯到Fredholm，Hilbert，Nystrom，Hadamard等大数学家 
* 在20世纪初，因为变分方法的崛起，积分方程方法在分析中的重要性减弱了，因为经典的积分方程的理论缺乏存在唯一性的精确结果
* 20世纪中叶，积分方程方法又随着数值PDE求解的需求而再次获得关注：
  * 相比于domain method(有限差分与有限元)，对于复杂的几何，简化了网格生成的难度。因为边界积分方程只需要生成边界的mesh，而不是volume的mesh
  * 无界域问题的优势
  * 对于参数依赖的高频问题（更稳定）
  * Better Condition Number compared to domain discretization
  * 难点：数值积分以及求解线性代数方程组--在80年代得到了高速的发展
  * 1980-1990 Galerkin方法被用以离散边界积分方程--之前用的是Nystrom以及配点法（相关的早期著作是Atkinson的关于第二类积分方程的书）
    * Galerkin方法的优势在于能够对一整大类(a very general class of)边界积分方程来分析稳定性、相容性以及收敛性，相关的著作有：
      * Hsiao & Wendland: Boundary Integral Eequations 2008
      * McLean: Strongly Elliptic Systems and Boundary Integral Equations 2000
      * Nedelec: Acoustic and Electromagnetic Equations 2001
    * 在实际应用中的Breakthrough，特别是对三维问题，来自于装配矩阵时积分的逼近以及边界积分算子（非局部）的表示的快速算法
* 问题一：关联椭圆边值问题的所有的边界积分方程是哪些？
* 问题二：第二类积分方程具有怎样的算子形式？

## 第一章：简例
* 静电场，方程由如下的假定给出：
  * 电荷分布为 rho，得到电位移的散度为 rho
  * 静电场无旋，因此 curl E = 0。E是有势场
  * 电位移与静电场差一个系数为介电常数
  * 假设电荷分布与一个有界区域 Omega 中，将区域分解成 Omega 的内部与外部，假定边界 Gamma 具有充分的光滑性以及介电常数光滑（使得散度定义合理），并且介电常数在外部为常数，则在 Omega 之外为一个Laplace方程（介电常数是scalar）
  * 边界条件的讨论：
    * 尤其是解在边界上的连续性条件
    * 在3.3节中，从数学的角度解释了（Check)，静电场E需要在界面处切向分量连续，而电位移D则需要在界面处法向分量连续（注意到D = epsilon E）
    * 结合potential approach，D的连续性可以通过potential flux = epsilon ∂Phi/∂n 的连续性得到
    * E的连续性即Phi沿着边界Gamma的导数连续，内外均可以有一个常数，选取某点使得在边界上内外的Phi相同，则由于切向导数均相同，得到Phi在界面内外的连续
    * **总结：势能在界面上连续，在法向上potential flux连续**
    * 如果epsilon是张量呢？
  * 无界域问题的无穷远衰减条件：势能以 r^-1 衰减，gradient以平方衰减
* 求解区域分解的Interface问题的积分方程方法
  * 求解内部方程：（rho紧支）利用全空间Laplace方程的基本解
    * 基本解： G(z) ~ 1/|z| (R3), 满足 ∆G= delta
    * 因此对于R3中的源项 rho，解为 Newton potential：N(x)=∫ G(x-y) rho(y)，形式上在 ∆ 作用下，变成delta函数和rho的卷积，从而满足方程
    * 在rho紧支的情况下，Newton potential满足衰减条件
    * 注：这样的解在无穷远满足decay条件，但是在interface上未必满足Dirichlet边界条件（尚且不论flux连续条件）
  * 但总之整体的解可以写成N加上另一函数（或许仅仅定义在外区域），显然这个函数在外区域上满足Laplace方程，修改过的Dirichelt边界条件（原来的减去N的）以及无穷远的衰减条件
    * 由于在内部没有source，所以待定边界上的密度sigma，将解表示成 Phi0(x) = ∫_Gamma G(x-y) sigma(y) (一个边界上的积分)
    * 上述定义的被称作single layer potential S(sigma)
    * 单层位势满足除了Interface以外的Laplace方程（这解释了在内部加上它不会影响内部的Poisson方程），同时对sigma在Gamma上连续满足无穷远衰减条件
    * 所以只要通过Phi0在边界上等于修改过的Dirichlet边界条件解出density sigma即可
    * 上述定义的单层位势定义在Interface之外，可以证明（3.1.16），它可以连续地延拓到界面处
    * **总结1：求解Dirichlet边界的Poisson方程，先用Newton Potential处理source，再用单层位势处理剩下的Dirichelt边界条件**
* G(x,y)可以求各阶偏导数，只要不要乘以x的函数，作为双变量（x，y）的函数仍然在x的∆作用下为0，由此可以定义双层位势
  * kernel： k(x,y) = <n(y), ∂yG(x-y)>
  * 对应于density：theta的双层位势：D(theta)(x) = ∫_Gamma k(x,y) theta(y), x out of Gamma
  * 连续延拓到Gamma上，记x取值在Gamma上为 K(theta)(x)，则延拓为 -1/2 theta(x) + K(theta)(x) （显然意味着D(theta)(x)在Gamma上并不连续），K是Gamma以外连续延拓到Gamma上的（单层位势没有这个问题）
  * 由于theta既出现在积分之内，又出现在积分之外，被称作第二类边界积分方程


## Chap4: 边界积分方程连续情形的存在唯一性迁移到Galerkin离散上


## Chap5: 计算细节：组装矩阵
* case by case quadrature method
  * singularity of the integrand: Duffy coordinates (Transformation) [83]

## Chap6: Iteration Method
* 经典的迭代方法的收敛由矩阵的条件数决定，三类情形：
  * 非对称的系统，但是有有界的条件数  —————— MinRes，迭代次数只与条件数相关
  * 对称正定的系统并且条件数以 √N 增长，N为边界有限元空间的维数，此外边界积分算子是smoothing的（作用后的光滑性比之前更高一阶） 
  * 边界积分算子不具有光滑性 （与椭圆边值问题的有限元离散紧密关联，使用多重网格）
  * 余下两种使用  ——————  cg method，迭代次数以 N 的1/4幂次增长
* 