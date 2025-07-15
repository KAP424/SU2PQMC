# 无序算符的构造

Ising model 的$Z_2$ 对称性：
$$
X_A=\Pi_{i\in A}\sigma^x_i
$$
repulsive Hubbard model 的U(1)对称性：每个site变化相同的相位不变，不同site给不同的相位等同于产生disorder
$$
X^c_{A}(\alpha)=\Pi_{i\in A} e^{i\alpha (n_{i\uparrow}+n_{i\downarrow})}
$$
attractive Hubbard model 的U(1)对称性
$$
X^\sigma_{A}(\alpha)=\Pi_{i\in A} e^{i\alpha (n_{i\uparrow}-n_{i\downarrow})}
$$




# 无序算符DQMC推导

$$
X_M(\alpha) = \frac{\sum \det[P^T B_s(2\theta,\theta) e^T B_s(\theta,0) P]}{\sum \det[P^T B_s(2\theta,0) P]}= \frac{\langle \Psi | U(2\theta,\theta) \prod_{i\in n} e^{i\alpha n_i} U(\theta,0) | \Psi \rangle}{\langle \Psi | U(2\theta,0) | \Psi \rangle}
$$

$$
= \sum_s \frac{\det[P^T B_s(2\theta,0) P]}{\sum \det[P^T B_s(2\theta,0) P]} \frac{\det[P^T B_s(2\theta,\theta) e^T B_s(\theta,0) P]}{\det[P^T B_s(2\theta,\theta) P]}= \sum_s P_s X_{M,s}^c
$$

$$
\quad (I - G) = BCA \quad B= B_s(\tau,0) P \quad C=[P^\dagger B_s(2\theta,0) P]^{-1} A=\quad P^\dagger B_s(2\theta,\tau)]
$$

e.g. 	A 5x10 B 10x5 C 5x5
$$
G_s(\tau) = [I - B_s(\tau,0) P [P^\dagger B_s(2\theta,0) P]^{-1} P^\dagger B_s(2\theta,\tau)]
$$

$$
X_{M,s} = \det \left\{ P^\dagger B_s(2\theta,\theta) e^T B_s(\theta,0) P [P^\dagger B_s(2\theta,0) P]^{-1} \right\}
$$


$$
e^T = \Delta + I \qquad det(I+AB) = det(I+BA)
$$

$$
\Delta = \begin{pmatrix} e^{i\alpha}-1 & & \\ & \ddots & \\ & & e^{i\alpha}-1 \\& & &0 \\ & & && \ddots\\ & & &&&0 \end{pmatrix}
$$

$$
X_{M,s}= \det\{I +\textcolor{red}{ P^\dagger B_s(2\theta,\theta)} \Delta B_s(\theta,0) P[P^\dagger B_s(2\theta,0) P]^{-1}\}
$$

$$
= \det\{I + \Delta B_s(0,0) P[P^\dagger B_s(2\theta,0) P]^{-1} \textcolor{red}{ P^\dagger B_s(2\theta,\theta)}\}
$$

$$
= \det\{I + \Delta [I - G_s(\theta)]\}= \det\{I + e^{i\alpha} [I - G_{M,s}(\theta)]\}\\
=det\{G_{M,s}+[I-G_{M,s}(\theta)]e^{i\alpha} \}
$$
# 增量version

$$
X_{M}=\frac{\sum_{s} P_{s} X_{M, s}} {\sum_{s} P_{s}}=\frac{\sum_{s} P_{s} ( X_{M,s} )^{\frac{1}{N}}} {\sum_{s} P_{s}} \cdots \frac{\sum_{s} P_{s} ( X_{M,s} )^{\lambda+\frac{1}{N}}} {\sum_{s} P_{s} ( X_{M,s} )^{\lambda}} \cdots \frac{\sum_{s} P_{s} X_{M,s}} {\sum_{s} P_{s} ( X_{M,s} )^{\lambda-\frac{1}{N}}}
$$
$$
\mathcal{P}_{s}^{\lambda}=P_{s} \left( X_{M,s}\right)^{\lambda}, \frac{P_{s^{\prime}}} {P_{s}}=1+\Delta_{x} ( 1-G_{xx} )=R
$$
$$
\mathcal{R}^\lambda=\frac{\mathcal{P}_{s^{\prime}}^{\lambda}} {\mathcal{P}_{s}^{\lambda}}=\frac{P_{s^{\prime}}} {P_{s}} \operatorname* {det} [ \frac{X_{M, s^\prime}} {X_{M, s}} ]^\lambda
$$
由于
$$
G_{s^{\prime}}(\theta)=G_{s}(\theta)+\frac{\Delta} {R} G_{s}( \theta, \tau )  [:, x] \otimes G_{s} ( \tau,\theta ) [:, s ]
$$
有
$$
\frac{X_{M,s^\prime}}  {X_{M,s}}=\frac{G_{M,s^\prime}+e^{i\alpha} ( I-G_{M,s\prime} )} {G_{M,s}+e^{i\alpha} ( I-G_{M,s} )}=\frac{e^{i\alpha}+( 1-e^{\alpha} ) G_{M,s^\prime}} {e^{\theta}+( 1-e^{-\theta} ) G_{M,s}}
$$
$$
= I+( 1-e^{i \theta} ) \frac{\Delta}{R} G ( \theta, \tau)[:,x] \otimes G (\tau,\theta) [x,:]X_{M,s}^{-1}
$$
$$
\vec{a}=G(\theta,\tau)[:,x] \qquad \vec{b}=G(\tau,\theta)[x,:]X_{M,s}^{-1} \qquad \Gamma=\vec{b}\cdot\vec{a}
$$

所以$(det[I+AB]=det[I+BA])$
$$
det[\frac{X_{M,s^\prime}}  {X_{M,s}}]=1+(1-e^{i\alpha})\frac{\Delta}{R}\Gamma
$$

$$
\mathcal{R}^\lambda=R\cdot [1+(1-e^{i\alpha})\frac{\Delta}{R}\Gamma]^\lambda
$$

再根据$[ A+u \otimes v ]^{-1}=A^{-1}-\frac{( A^{-1} u ) \otimes( v A^{-1} )} {1+v \cdot A^{-1} u} $有
$$
[\frac{X_{M,s^\prime}}  {X_{M,s}}]^{-1}=I-\frac{(1-e^{i\alpha})\Delta}{[R+(1-e^{i\alpha})\Delta]\Gamma} \vec{a}\otimes\vec{b}
$$


------------

Particle-Hole Symmetry下
$$
G^\uparrow=<\vec{c_\uparrow}\vec{c_\uparrow}^\dagger>
\\
G^\uparrow=<\vec{h^\dagger}\vec{h}>=1-<\vec{h}\vec{h^\dagger}>=1-G^{\uparrow*}
$$

$$
X^\uparrow(\alpha)=G^\uparrow+e^{i\alpha}(I-G^\uparrow)
\\
X^\downarrow(\alpha)=G^\downarrow+e^{i\alpha}(I-G^\downarrow)=(I-G^{\uparrow*})+e^{i\alpha}G^{\uparrow*}=e^{i\alpha}[G^{\uparrow*}+e^{-i\alpha}(I-G^{\uparrow*})]
\\
\Rightarrow X^\downarrow(\alpha)=e^{i\alpha} X^{\uparrow*}(\alpha)
$$

故而charger 无序算符相位固定，计算无符号问题：
$$
X^c_{M,s}(\alpha)=X^\uparrow_{M,s}(\alpha)X^\downarrow_{M,s}(\alpha)=e^{i\alpha}|X^\uparrow_{M,s}(\alpha)|^2 \\
$$

> 其他对称性下的无序算符怎么构造
>
> 及其对应的无符号问题的DQMC计算



