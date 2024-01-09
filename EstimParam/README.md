Given
- T observations $Z_1, \cdots, Z_T$ taking values in $\mathbb{Z}_{\geq 0}$
- T mean values $\Phi_1, \cdots, \Phi_T$ taking values on $\mathbb{R}_{\geq 0}$
- Two initial values for the reproduction number $R_1, R_2$, taking values in $\mathbb{R}_{>0}$
- Two parameters $\lambda_R$ and $\lambda_O$ taking values in $\mathbb{R}_{>0} \times \mathbb{R}_{>0}$

two target distributions are defined.

**Model "no mixture"**
A target distribution $\pi$ on $(\mathbb{R}_{>0} \times \mathbb{R})^{T-2}$. The log-density $\ln \pi$ is given by

$$ 
\begin{split}
(R_3, O_3, \cdots, R_T, O_T) & \mapsto \sum_{t=3}^T \Bigl( Z_t \ln( R_t \Phi_t+O_t) - (R_t \Phi_t + O_t) \Bigr)  \\
& - \frac{\lambda_R}{4} \sum_{t=3}^T | R_t - 2 R_{t-1} + R_{t-2} | + (T-2) \ln \lambda_R   \\
& - \lambda_0  \sum_{t=3}^T |O_t| + (T-2) \ln \lambda_0.
\end{split}
$$

up to an additive constant. The support of the distribution is for $t=3, \cdots, T$: 

$$
R_t > 0; \qquad \qquad  R_t \Phi_t + O_t > 0  \quad \text{if} \quad Z_t >0, \quad \text{or} \quad R_t \Phi_t + O_t \geq 0  \quad \text{if} \quad Z_t  \geq 0.
$$

**Model "mixture"**

