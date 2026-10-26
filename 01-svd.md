# 作业一：奇异值分解（SVD）的过程与实例

## 1. 什么是奇异值分解

**奇异值分解（Singular Value Decomposition, SVD）** 是线性代数中最重要的矩阵分解之一。对任意矩阵 $A \in \mathbb{R}^{m \times n}$（不要求方阵、不要求可逆、不要求对称），都存在分解：

$$
A = U\Sigma V^{T}
$$

其中：

- $U \in \mathbb{R}^{m \times m}$ 是**正交矩阵**（$U^{T}U = I$），其列向量 $u_1,\dots,u_m$ 称为 $A$ 的**左奇异向量**，是 $AA^{T}$ 的单位特征向量；
- $\Sigma \in \mathbb{R}^{m \times n}$ 是**对角矩阵**，对角元素 $\sigma_1 \ge \sigma_2 \ge \cdots \ge \sigma_r > 0$ 称为**奇异值**，$r = \mathrm{rank}(A)$，其余位置全为 0；奇异值是 $A^{T}A$（同时也是 $AA^{T}$）特征值的算术平方根；
- $V \in \mathbb{R}^{n \times n}$ 是**正交矩阵**，其列向量 $v_1,\dots,v_n$ 称为**右奇异向量**，是 $A^{T}A$ 的单位特征向量。

> **几何意义**：SVD 告诉我们，任意线性变换 $x \mapsto Ax$ 都可以分解为三步——先在输入空间旋转/反射（$V^T$），再沿各坐标轴缩放 $\sigma_i$（$\Sigma$，会把 $n$ 维空间"压扁"到 $r$ 维），最后在输出空间旋转/反射（$U$）。这把一个"任意扭曲"的变换拆解成了刚体运动 + 拉伸。

## 2. SVD 的存在性原理（推导思路）

为什么这样的分解一定存在？核心思路如下：

1. $A^{T}A$ 是**半正定对称矩阵**，故可正交对角化：存在正交阵 $V$ 使  
   $V^{T}(A^{T}A)V = \mathrm{diag}(\lambda_1,\dots,\lambda_n),\quad \lambda_1\ge\cdots\ge\lambda_r>0=\lambda_{r+1}=\cdots$
2. 取 $\sigma_i = \sqrt{\lambda_i}$，记 $V = [v_1,\dots,v_n]$。对 $i \le r$，定义  
   $u_i = \frac{1}{\sigma_i}Av_i$
3. 验证 $\{u_1,\dots,u_r\}$ 是标准正交的：  
   $u_i^{T}u_j = \frac{v_i^{T}A^{T}Av_j}{\sigma_i\sigma_j} = \frac{\lambda_j}{\sigma_i\sigma_j}v_i^{T}v_j = \delta_{ij}$
4. 将其扩充为 $\mathbb{R}^m$ 的标准正交基，得到正交阵 $U$，并直接验证 $Av_i = \sigma_i u_i$ 对所有 $i$ 成立（$i>r$ 时两边都是 0），于是 $AV = U\Sigma$，即 $A = U\Sigma V^{T}$。∎

## 3. 手动求 SVD 的标准流程

对 $A \in \mathbb{R}^{m\times n}$（设 $m \ge n$），按以下五步操作：

**Step 1** 计算 Gram 矩阵 $A^{T}A$（$n \times n$ 对称半正定）。

**Step 2** 求 $A^{T}A$ 的特征值，从大到小排序：$\lambda_1 \ge \cdots \ge \lambda_n \ge 0$，奇异值 $\sigma_i = \sqrt{\lambda_i}$，填入 $\Sigma$ 的前 $r$ 个对角元（$r$ 为非零特征值个数）。

**Step 3** 对每个 $\lambda_i$ 求单位特征向量 $v_i$，组成 $V = [v_1,\dots,v_n]$。

**Step 4** 对 $i = 1,\dots,r$ 计算左奇异向量 $u_i = \dfrac{Av_i}{\sigma_i}$；若 $r < m$，用 Gram–Schmidt 补全 $u_{r+1},\dots,u_m$（它们是 $AA^{T}$ 对应零特征值的特征向量，即 $A$ 列空间的正交补）。

**Step 5** 拼出 $U, \Sigma, V^{T}$，并用 $A \stackrel{?}{=} U\Sigma V^{T}$ 验证。

## 4. 手算实例

取矩阵

$$
A = \begin{bmatrix} 1 & 0 \\ 0 & 1 \\ 1 & 1 \end{bmatrix} \in \mathbb{R}^{3\times 2}
$$

下面严格按五步流程分解。

### Step 1：计算 $A^{T}A$

$$
A^{T}A = \begin{bmatrix} 1 & 0 & 1 \\ 0 & 1 & 1 \end{bmatrix}\begin{bmatrix} 1 & 0 \\ 0 & 1 \\ 1 & 1 \end{bmatrix} = \begin{bmatrix} 2 & 1 \\ 1 & 2 \end{bmatrix}
$$

### Step 2：求特征值与奇异值

特征方程：

$$
\det(A^{T}A - \lambda I) = \begin{vmatrix} 2-\lambda & 1 \\ 1 & 2-\lambda \end{vmatrix} = (2-\lambda)^2 - 1 = \lambda^2 - 4\lambda + 3 = 0
$$

解得 $\lambda_1 = 3,\ \lambda_2 = 1$（均大于 0，故 $\mathrm{rank}(A)=2$）。奇异值：

$$
\sigma_1 = \sqrt{3}, \qquad \sigma_2 = 1
$$

于是

$$
\Sigma = \begin{bmatrix} \sqrt{3} & 0 \\ 0 & 1 \\ 0 & 0 \end{bmatrix}_{3\times 2}
$$

### Step 3：求右奇异向量 $V$

- 对 $\lambda_1 = 3$：$(A^{T}A - 3I)x = \begin{bmatrix} -1 & 1 \\ 1 & -1 \end{bmatrix}x = 0 \Rightarrow x_1 = x_2$，单位化得  
  $v_1 = \frac{1}{\sqrt{2}}\begin{bmatrix} 1 \\ 1 \end{bmatrix}$
- 对 $\lambda_2 = 1$：$(A^{T}A - I)x = \begin{bmatrix} 1 & 1 \\ 1 & 1 \end{bmatrix}x = 0 \Rightarrow x_1 = -x_2$，单位化得  
  $v_2 = \frac{1}{\sqrt{2}}\begin{bmatrix} 1 \\ -1 \end{bmatrix}$

（可验证 $v_1 \perp v_2$，$V$ 为正交阵。）

$$
V = \frac{1}{\sqrt{2}}\begin{bmatrix} 1 & 1 \\ 1 & -1 \end{bmatrix}, \qquad V^{T} = \frac{1}{\sqrt{2}}\begin{bmatrix} 1 & 1 \\ 1 & -1 \end{bmatrix}
$$

### Step 4：求左奇异向量 $U$

$$
u_1 = \frac{Av_1}{\sigma_1} = \frac{1}{\sqrt{3}}\cdot\frac{1}{\sqrt{2}}\begin{bmatrix} 1 \\ 0 \\ 1 \end{bmatrix} + \frac{1}{\sqrt{3}}\cdot\frac{1}{\sqrt{2}}\begin{bmatrix} 0 \\ 1 \\ 1 \end{bmatrix} \ \Rightarrow\ u_1 = \frac{1}{\sqrt{6}}\begin{bmatrix} 1 \\ 1 \\ 2 \end{bmatrix}
$$

$$
u_2 = \frac{Av_2}{\sigma_2} = \frac{1}{\sqrt{2}}\begin{bmatrix} 1 \\ 0 \\ 1 \end{bmatrix} - \frac{1}{\sqrt{2}}\begin{bmatrix} 0 \\ 1 \\ 1 \end{bmatrix} \ \Rightarrow\ u_2 = \frac{1}{\sqrt{2}}\begin{bmatrix} 1 \\ -1 \\ 0 \end{bmatrix}
$$

因 $r = 2 < m = 3$，还需补一个与 $u_1, u_2$ 都正交的单位向量。取叉积再单位化：

$$
u_1 \times u_2 \propto \begin{vmatrix} i & j & k \\ 1 & 1 & 2 \\ 1 & -1 & 0 \end{vmatrix} = (2,\ 2,\ -2) \ \Rightarrow\ u_3 = \frac{1}{\sqrt{3}}\begin{bmatrix} 1 \\ 1 \\ -1 \end{bmatrix}
$$

于是

$$
U = \begin{bmatrix} \dfrac{1}{\sqrt{6}} & \dfrac{1}{\sqrt{2}} & \dfrac{1}{\sqrt{3}} \\[6pt] \dfrac{1}{\sqrt{6}} & -\dfrac{1}{\sqrt{2}} & \dfrac{1}{\sqrt{3}} \\[6pt] \dfrac{2}{\sqrt{6}} & 0 & -\dfrac{1}{\sqrt{3}} \end{bmatrix}
$$

### Step 5：验证 $A = U\Sigma V^{T}$

先算 $\Sigma V^{T}$：

$$
\Sigma V^{T} = \begin{bmatrix} \sqrt{3} & 0 \\ 0 & 1 \\ 0 & 0 \end{bmatrix}\cdot\frac{1}{\sqrt{2}}\begin{bmatrix} 1 & 1 \\ 1 & -1 \end{bmatrix} = \frac{1}{\sqrt{2}}\begin{bmatrix} \sqrt{3} & \sqrt{3} \\ 1 & -1 \\ 0 & 0 \end{bmatrix}
$$

再左乘 $U$，第 1 列：

$$
U\Sigma V^{T}\text{ 的第1列} = \frac{\sqrt{3}}{\sqrt{2}}u_1 + \frac{1}{\sqrt{2}}u_2 = \frac{\sqrt{3}}{\sqrt{2}}\cdot\frac{1}{\sqrt{6}}\begin{bmatrix}1\\1\\2\end{bmatrix} + \frac{1}{\sqrt{2}}\cdot\frac{1}{\sqrt{2}}\begin{bmatrix}1\\-1\\0\end{bmatrix} = \frac{1}{2}\begin{bmatrix}1\\1\\2\end{bmatrix} + \frac{1}{2}\begin{bmatrix}1\\-1\\0\end{bmatrix} = \begin{bmatrix}1\\0\\1\end{bmatrix} \checkmark
$$

第 2 列同理可得 $\begin{bmatrix}0\\1\\1\end{bmatrix}$，正是 $A$ 的第 2 列。验证通过。

> **数值验证**：用 NumPy 的 `numpy.linalg.svd` 计算，奇异值为 $1.732051 = \sqrt{3}$ 和 $1$，与手算完全一致（$U, V$ 的符号可能整体取反，不影响分解）。

### 4.1 奇异值的另一种验证：Frobenius 范数守恒

Frobenius 范数的平方等于奇异值平方和：

$$
\|A\|_F^2 = 1^2+0^2+0^2+1^2+1^2+1^2 = 4 = \sigma_1^2 + \sigma_2^2 = 3 + 1 \checkmark
$$

## 5. SVD 的重要性质与应用：最优低秩逼近

SVD 可写成求和形式 $A = \sum_{i=1}^{r} \sigma_i u_i v_i^{T}$（每项 $\sigma_i u_i v_i^T$ 都是秩 1 矩阵）。**截断 SVD** 只保留前 $k$ 项：

$$
A_k = \sum_{i=1}^{k} \sigma_i u_i v_i^{T}
$$

由 **Eckart–Young 定理**，$A_k$ 是所有秩不超过 $k$ 的矩阵中对 $A$ 的最优逼近（Frobenius 范数意义下），且误差

$$
\|A - A_k\|_F^2 = \sum_{i=k+1}^{r}\sigma_i^2
$$

**在本例中取 $k=1$**：

$$
A_1 = \sigma_1 u_1 v_1^{T} = \sqrt{3}\cdot\frac{1}{\sqrt{6}}\begin{bmatrix}1\\1\\2\end{bmatrix}\cdot\frac{1}{\sqrt{2}}\begin{bmatrix}1 & 1\end{bmatrix} = \begin{bmatrix} 0.5 & 0.5 \\ 0.5 & 0.5 \\ 1 & 1 \end{bmatrix}
$$

$$
\|A - A_1\|_F^2 = \sigma_2^2 = 1
$$

即只用 **3 个奇异值分量**（$A_1$ 的自由度为 $1\times(3+2)=5$ 个参数）就以 $\frac{3}{4}$ 的能量保留率逼近了原矩阵（能量占比 $\frac{\sigma_1^2}{\sigma_1^2+\sigma_2^2} = \frac{3}{4}$）。

**应用举例**：

| 应用         | 原理                                       |
| ---------- | ---------------------------------------- |
| 图像压缩 / PCA | 保留前 $k$ 个大奇异值，丢弃小奇异值对应的噪声分量              |
| 数据降维       | 用 $V$ 的前 $k$ 列做投影（等价于 PCA）               |
| 伪逆与最小二乘    | $A^{+} = V\Sigma^{+}U^{T}$，病态/欠定方程组的稳定求解 |
| 推荐系统       | 矩阵补全（如 Funk-SVD）把评分矩阵分解为用户/物品隐向量         |
| 语义分析       | LSA 对词-文档矩阵做截断 SVD 提取主题                  |

## 6. 小结

- SVD 对**任意矩阵**存在且唯一（奇异值确定，$U,V$ 相差符号/退化子空间内的正交变换）。
- 计算套路固定：$A^{T}A \to$ 特征值 $\to$ 奇异值 $\to V \to u_i = Av_i/\sigma_i \to$ 补全 $U \to$ 验证。
- 奇异值刻画了矩阵的"能量分布"，截断 SVD 给出最优低秩逼近，是压缩、降噪、降维的统一工具。
