$a>0$

-   二分量：

  - 排斥U: $a=\Delta U$
    $$
    e^{-aO^2}=\frac{1}{2} \sum_{s=\pm 1} e^{i\alpha sO}
    \\
    \cos\alpha = e^{-a}\qquad
    $$
  
- 吸引U: $a=-\Delta U$

$$
e^{aO^2}=\frac{1}{2} \sum_{s=\pm 1} e^{\alpha sO}
\\
\cosh\alpha = e^{a}
$$

- 四分量

$$
e^{-aO^2}=
\frac{1}{4} \sum_{\substack{l=\pm 1 , \pm 2}} \gamma(l) e^{i\sqrt{a} \eta(l) O} 
\\
e^{aO^2}=
\frac{1}{4} \sum_{\substack{l=\pm 1 , \pm 2}} \gamma(l) e^{-\sqrt{a} \eta(l) O}\\
\gamma(\pm 1) = \left(1+\frac{\sqrt{6}}{3}\right) \\
\gamma(\pm 2) = \left(1-\frac{\sqrt{6}}{3}\right) \\

\eta(\pm 1) = \pm\sqrt{2(3-\sqrt{6})} \\
\eta(\pm 2) = \pm\sqrt{2(3+\sqrt{6})} \\
$$

```julia
using LinearAlgebra
using Plots

a=1


α=acos(exp(-a))
C=1/2

x=collect(-2:0.01:2)
dy=zeros(length(x))
for i in eachindex(x)

    dy[i]=exp(-a*x[i]^2)-C*exp(1im*α*x[i])-C*exp(-1im*α*x[i])

end
scatter(x,dy,label="Two component ")

savefig("1.png")


a=1
α=sqrt(a)
γ=[1+sqrt(6)/3 1+sqrt(6)/3 1-sqrt(6)/3 1-sqrt(6)/3]
η=[sqrt(2*(3-sqrt(6))) -sqrt(2*(3-sqrt(6))) sqrt(2*(3+sqrt(6))) -sqrt(2*(3+sqrt(6)))]
C=1/4

x=collect(-2:0.01:2)
dy=zeros(length(x))*1im

for i in eachindex(x)
    dy[i]=exp(-a*x[i]^2)
    for j in 1:4
        dy[i]-=C*γ[j]*exp(1im*sqrt(a)*η[j]*x[i] )
    end
end

scatter(x,abs.(dy),label="Four component ")
savefig("2.png")

```

<img src="E:\桌面\JuliaDQMC\basic\HS\1.png" alt="1" style="zoom:50%;" />

<img src="E:\桌面\JuliaDQMC\basic\HS\2.png" alt="2" style="zoom:50%;" />

排斥的有两个channel

![image-20250630114336221](C:\Users\KAP\AppData\Roaming\Typora\typora-user-images\image-20250630114336221.png)





# 半满填充的hubbard 模型

## repulsive-U

$$
H=t \sum_{\langle i, j \rangle} ( c_{i, \sigma}^{\dagger} c_{j, \sigma}+h. c. )+U\sum_i n_{i\uparrow}n_{i\downarrow}\\
H_{int}=U\sum_i n_{i\uparrow}n_{i\downarrow}=\frac{U}{2}\sum_i(n_{i\uparrow}+n_{i\downarrow}-1)^2
$$

虚实切片
$$
e^{-\Delta H_{int}}=\Pi e^{-\frac{\Delta U}{2} (n_{i\uparrow}+n_{i\downarrow}-1)^2}=\frac{1}{4}\sum_{\vec{s}}\Pi_i \gamma(s_i)e^{i\sqrt{\frac{\Delta U}{2} }\eta(s_i)[n_{i\uparrow}+n_{i\downarrow}-1] }
$$
Particle-Hole 变换：$\cases{c_{i\uparrow} \rightarrow c_i\\ c_{i\downarrow} \rightarrow (-1)^i h_i^\dagger}$
$$
e^{-\Delta H_{int}}=\frac{1}{4}\sum_{\vec{s}}\Pi_i \gamma(s_i)e^{i\sqrt{\frac{\Delta U}{2} }\eta(s_i)[n_{i\uparrow}-h_i^\dagger h_i] }
$$
从而有：$<c_{i\uparrow}c_{j\uparrow}^\dagger>=<h_ih_j^\dagger>^\star$，又$G_\uparrow=[<c_{i\uparrow}c_{j\uparrow}^\dagger>]$
$$
G_\downarrow=[<c_{i\downarrow}c_{j\downarrow}^\dagger> ]=[(-1)^{i+j}<h_i^\dagger h_j>]=I-[(-1)^{i+j}<h_{j}h_{i}^\dagger>]
\\
=I-[(-1)^{i+j}<c_{j\uparrow}c_{i\uparrow}^\dagger>^\star]=I-[(-1)^{i+j}G_\uparrow^\dagger]
$$

由$[(-1)^{i+j}A]=U\times A\times U \qquad U=\begin{bmatrix}
1 &  & \\
 & -1 & \\ & & -1\\&&&\ddots
\end{bmatrix}$,有：$G_\downarrow=I-UG_\uparrow^\dagger U$

> 后文的G都是$ G(\theta) $ 的缩写

1. 纠缠熵：
   $$
   e^{-S_A}=\sum_{s_1s_2} P_{s_1s_2} \det g_{s_1s_2}^A
   \\
   g_{s_1s_2}^A=g_{s_1s_2\uparrow}^A\otimes g_{s_1s_2\downarrow}^A 
   \\
   g_{s_1s_2\uparrow}^A=G_{s_1\uparrow}^AG^A_{s_2\uparrow}+(I-G^A_{s_1\uparrow})(I-G^A_{s_2\uparrow}) 
   \\
   g_{s_1s_2\downarrow}^A=G^A_{s_1\downarrow}G^A_{s_2\downarrow}+(I-G^A_{s_1\downarrow})(I-G^A_{s_2\downarrow})
   $$
   

   $$
   g_{s_1s_2\downarrow}^A=U [{G^A_{s_1\uparrow}}^\dagger {G^A_{s_2\uparrow}}^\dagger+(I-{G^A_{s_1\uparrow}}^\dagger)(I-{G^A_{s_2\uparrow}}^\dagger)] U
   
   =U{g_{s_2s_1\uparrow}^A}^\dagger U
   $$
    从而：
   $$
   \det g_{s_1s_2\downarrow}^A = \det {g_{s_2s_1\uparrow}^A}^\dagger=(\det g_{s_2s_1\uparrow}^A)^\star \neq (\det g_{s_1s_2\uparrow}^A)^\star
   $$
   <u>**BUT，数值结果表明：**</u>
   $$
   \det g_{s_2s_1\uparrow}^A = \det g_{s_1s_2\uparrow}^A\\
   $$

   > 但同时 $  g_{s_2s_1\uparrow}^A \neq g_{s_1s_2\uparrow}^A , ({g_{s_1s_2\uparrow}^A})^T$  $G^A_{s_1\uparrow}G^A_{s_2\uparrow}\neq G^A_{s_2\uparrow}G^A_{s_1\uparrow}$
   >
   > ？？？？？？
   >
   > 这个在数值上为什么会相等也十分奇怪

   从而：
   $$
   \det g_{s_1s_2\downarrow}^A=(\det g_{s_1s_2\uparrow}^A)^\star
   \\
   \det g_{s_1s_2}^A= |\det g_{s_1s_2\uparrow}^A|^2
   $$

2. charge disorder operator

$$
X_{A,s}^{c}(\alpha)=X_{A,s\uparrow}(\alpha)\otimes X_{A,s\downarrow}(\alpha)
\\
X_{A,s\uparrow}(\alpha)=\det\{G^A_{s\uparrow}+e^{i\alpha}[I-G^A_{s\uparrow}] \}
\\
X_{A,s\downarrow}(\alpha)=\det\{G^A_{s\downarrow}+e^{i\alpha}[I-G^A_{s\downarrow}] \}
=\det\{I-{G^A_{s\uparrow}}^\dagger+e^{i\alpha}{G^A_{s\uparrow}}^\dagger \}
\\
=e^{i \alpha N} \det \{ {G^A_{s\uparrow}}^\dagger+e^{-i\alpha}[I-{G^A_{s\uparrow}}^\dagger] \}=e^{i\alpha N} {X_{A,s\uparrow}(\alpha)}^\star
$$

从而：
$$
|X_{A,s}^{c}(\alpha)|=|X_{A,s\uparrow}(\alpha)|^2
$$


## attractive-U

$$
H=t \sum_{\langle i, j \rangle} ( c_{i, \sigma}^{\dagger} c_{j, \sigma}+h. c. )-U\sum_i n_{i\uparrow}n_{i\downarrow}
\\
H_{int}=-U\sum_i n_{i\uparrow}n_{i\downarrow}=\frac{U}{2}\sum_i[(n_{i\uparrow}-n_{i\downarrow})^2-1]\propto  \frac{U}{2}\sum_i(n_{i\uparrow}-n_{i\downarrow})^2
$$

虚实切片
$$
e^{-\Delta H_{int}}=\Pi e^{-\frac{\Delta U}{2} (n_{i\uparrow}-n_{i\downarrow})^2}=\frac{1}{4}\sum_{\vec{s}}\Pi_i \gamma(s_i)e^{i\sqrt{\frac{\Delta U}{2} }\eta(s_i)[n_{i\uparrow}-n_{i\downarrow}] }
$$
有：$G_\downarrow=G_\uparrow^\star$，从而：

1. 纠缠熵
   $$
   e^{-S_A}=\sum_{s_1s_2} P_{s_1s_2} \det g_{s_1s_2}^A
   \\
   g_{s_1s_2}^A=g_{s_1s_2\uparrow}^A\otimes g_{s_1s_2\downarrow}^A 
   \\
   g_{s_1s_2\downarrow}^A ={g_{s_1s_2\uparrow}^A}^\star
   $$

   从而
   $$
   \det g_{s_1s_2}^A= |\det g_{s_1s_2\uparrow}^A|^2
   $$

2. charge disorder operator
   $$
   X_{A,s\downarrow}(\alpha)=\det\{G^A_{s\downarrow}+e^{i\alpha}[I-G^A_{s\downarrow}] \}
   =\det\{{G^A_{s\uparrow}}^\star+e^{i\alpha}[I-{G^A_{s\uparrow}}^\star] \}
   \\
   ={X_{A,s\uparrow}(-\alpha)}^\star
   $$
   

   ​		从而:
   $$
   |X_{A,s}^{c}(\alpha)|=X_{A,s\uparrow}(\alpha)X_{A,s\downarrow}(\alpha)=X_{A,s\uparrow}(\alpha) {X_{A,s\uparrow}(-\alpha)}^\star
   $$

3. spin disorder operator
   $$
   X_{A,s}^{\sigma}(\alpha)=X_{A,s\uparrow}^{\sigma}(\alpha)X_{A,s\downarrow}^{\sigma}(-\alpha)=|X_{A,s\uparrow}^{\sigma}(\alpha)|^2
   $$
   


-----------

## Strange Thing

上面这些推导和数值似乎表明：排斥和吸引Hubbard模型的纠缠熵是对称的，即：$S_A(-U)=S_A(+U)$

如果纠缠熵可以表征相变，那排斥和吸引Hubbard模型的相图是对称的！？

但实际上排斥U对应的是反铁磁序，吸引U对应的是pairing序
