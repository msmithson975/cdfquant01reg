# Distribution Family

## Finite-Tailed CDF-Quantile Distributions

The CDF of the CDF-Quantile family (Smithson and Shou, 2017) is
``` math
G(x,\mu,\sigma) = F[U(H^{ - 1}(x),\mu,\sigma)],
```
for $`x \in \left[ {0,1} \right]`$, real location parameter $`\mu`$, and
dispersion parameter $`\sigma > 0`$. $`F`$ and $`H`$ are standard
symmetrical CDFs with support on the real line. the default $`U`$
function is
``` math
U(y,\mu ,\sigma) = (y - \mu)/\sigma
```
The finite-tailed subfamily comprises distributions that have defined
finite densities at 0 and at 1. These include any combination of the
standard Cauchy and Arcsinh distributions for $`F`$ and $`H`$, and the t
distribution with with 2 degrees of freedom for both $`F`$ and $`H`$
(i.e., Arcsinh-Arcsinh, Arcsinh-Cauchy, Cauchy-Arcsinh, Cauchy-Cauchy,
and t2-t2).

The finite-tailed subfamily, FTCDFQ (Smithson and Shou, 2024), employs
functions of a skew parameter, $`\theta`$. There are four alternative
locations for these functions, two with mappings
$`V:\left[ {0,1} \right] \to \left[ {0,1} \right]`$ applied to $`x`$ or
to $`F`$, and two with mappings $`W:\mathbb{R} \to \mathbb{R}`$ applied
to $`H^{-1}`$ or to $`U`$. Each pair’s function thus has an “inner” and
an “outer” position. The inner-position $`V`$ is applied to $`x`$ and
the outer-position $`V`$ is applied to $`F`$. Likewise, the
inner-position $`W`$ is applied to $`H^{-1}`$ and the outer-position
$`W`$ is applied to $`U`$. Given that the skew parameter $`\theta`$ also
influences location, 2-parameter as well as 3-parameter distributions
are viable, i.e., $`G\left( {x,\theta ,\sigma } \right)`$ and
$`G\left( {x,\mu,\theta ,\sigma } \right)`$.

$`V`$ is the so-called \`\`tilt’’ parameter function:
``` math
V(y,\theta) = e^{\theta}y/(e^{\theta}y + 1 - y),
```
for real $`\theta`$. It satisfies the requirements of being monotone
nondecreasing in $`y`$ and $`V(0,\theta) = 0`$ and $`V(1,\theta) = 1`$.
Its derivative has limits $`e^\theta`$ as $`y`$ goes to 0 and
$`e^{-\theta}`$ as $`y`$ goes to 1. $`W`$ is is a modified version of a
2-parameter transformation, $`\sinh( \theta  + \delta\sinh^{ - 1}(y))`$,
for real $`\theta`$ and $`\delta>0`$, where $`\theta`$ acts as a skew
parameter and $`\delta`$ controls tail thickness. Fixing $`\delta=1`$
gives
``` math
W(y,\theta) = \sinh(\theta  + \sinh^{ - 1}(y)),
```
whose derivative has limits $`e^\theta`$ as $`y`$ goes to $`\infty`$ and
$`e^{-\theta}`$ as $`y`$ goes to $`-\infty`$.

Smithson and Shou (2024) show that the inner-position $`W`$ and $`V`$
functions are equivalent when the $`H`$ function is ArcSinh, and
likewise the outer-position $`W`$ and $`V`$ functions are equivalent
when the $`F`$ function is ArcSinh. Thus, the FTCDFQ family has two sets
of 16 distinct distributions.

## Extended-Support Distributions

Kosmidis and Zeileis (2025) developed a method of extending the support
beyond 0 and 1 for distributions with support on the unit interval.
Their approach is applicable to any distribution whose support is
$`(0,1)`$ or $`[0,1]`$. They applied their method to the beta
distribution and implemented it in the betareg package. In this package
the method is applied to the 2-parameter FTCDFQ distributions. Let $`f`$
denote a pdf of $`x`$ whose support is the compact unit interval, with a
vector of parameters \$\bf{\beta}\$. Let $`F`$ denote the cdf, and
$`F^{-1}`$ the quantile function. An extended support pdf whose
distribution follows $`f`$ will have support on $`[-u, u+1]`$ for
$`u > 0`$. Denoting the extended support pdf by $`t`$, we define it as a
censored version of $`f`$ on $`[-u, u+1]`$:
``` math
\begin{array}{cl}
t(x,u,\beta) = F(u/(2u+1),\beta), x = 0\\
t(x,u,\beta) = f((u+x)/(2u+1),\beta),x \in \left( {0,1} \right)\\
t(x,u,\beta) = 1 - F((u+1)/(2u+1),\beta),x = 1
\end{array}
```
We also have the following identity for the extended-support
distribution $`t`$: \$\${T^{ - 1}}\left( {\gamma ,u,{\bf{\beta }}}
\right) = \left( {2u + 1} \right){F^{ - 1}}\left( {\gamma ,{\bf{\beta
}}} \right) - u\$\$ for $`\gamma  \in \left( {0,1} \right)`$, where
$`T^{-1}`$ denotes the quantile function for $`t`$. This equation
identifies the linear shift between the quantiles of $`f`$ and $`t`$. It
amounts to a principled way of either extending the boundaries from
$`[0,1]`$ to $`[-u, u +1]`$ or shrinking them to
$`[u/(2u+1), (u+1)/(2u+1)]`$.Thus, we can choose between an extended-
versus a shrunken-support interpretation of this censored-distribution
model.

That said, it is not always clear when we should prefer one over the
other. The extended-support interpretation infers from the observed
spikes at 0 and at 1 that we have underestimated the breadth of the
construct. One shrunken-support interpretation might be that people
whose responses are near the boundaries have “rounded off” their
responses to sit on the boundaries (e.g., someone whose subjective
probability assignment is 0.95 has rounded it off to 1). Another
interpretation is to infer that the scale is too sensitive to midrange
values at the sacrifice of sensitivity to extreme values and needs to be
compressed to make room for its genuine tails near 0 and 1.
