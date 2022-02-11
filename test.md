# Boundary Element Methods By Sauter and Schwab
* 积分方程方法具有悠久的历史，可以追溯到Fredholm，Hilbert，Nystrom，Hadamard等大数学家 
* 在20世纪初，因为变分方法的崛起，积分方程方法在分析中的重要性减弱了，因为经典的积分方程的理论缺乏存在唯一性的精确结果
* 20世纪中叶，积分方程方法又随着数值PDE求解的需求而再次获得关注：
    * 相比于domain method(有限差分与有限元)，对于复杂的几何，简化了网格生成的难度。因为边界积分方程只需要生成边界的mesh，而不是volume的mesh
    * 无界域问题的优势
    * 对于参数依赖的高频问题（更稳定）