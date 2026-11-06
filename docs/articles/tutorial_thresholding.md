# The effect of thresholding on correlations

## Theorem: The effect of thresholding on correlations between multivariate points

Given a finite set of m-dimensional real-valued vectors
S\subset\mathbb{R}^m, suppose:

- **S1.** There are n vectors in S and these vectors are indexed such
  that S can be treated as an n\times m matrix with the vectors s\in S
  as rows.
- **S2.** The columns s_j=\langle s\_{ij}\rangle\_{1\leq i\leq n} of S
  follow Gaussian distributions, i.e.,
  s\_{ij}\sim\mathcal{N}(\mu_j,\sigma_j).
- **S3.** The columns of S have pairwise correlations given by R, an
  m\times m symmetric matrix with ones along the diagonal, i.e.,
  R\_{j\_{1}j\_{2}}=\mathrm{corr}(s\_{j_1},s\_{j_2}).

And define:

- **D1.** The real-valued vectors \mu_S,\sigma_S\in\mathbb{R}^m hold,
  respectively, the column means \mu_j and column standard deviations
  \sigma_j of S, i.e., \mu_S=\langle\mu_j\rangle\_{1\leq j\leq m} and
  \sigma_S=\langle\sigma_j\rangle\_{1\leq j\leq m}.
- **D2.** The diagonal matrix D is the m\times m matrix with all zeros,
  except the diagonal, which contains \sigma_S, i.e., D\_{jj} = \sigma_j
  for all 1\leq j\leq m.
- **D3.** The matrix S^{\*} is the result of thresholding S into binary
  values of 0 and 1 according to the rule s^{\*}\_{ij}=1 if s\_{ij} \> 0
  and s^{\*}\_{ij}=0 otherwise.
- **D4.** For a square matrix M and integer vector v shorter than the
  sides of M, M\_{v} is the square matrix formed from the rows and
  columns indexed by v.

It follows that:

- **T1.** The covariance matrix for the columns of S is \Sigma = DRD.
- **T2.** The column means \mu^{\*}\_j of S^{\*} are given by the
  upper-tail Gaussian cumulative distribution function \Phi^{+}:
  \mu^{\*}\_j = \Phi^{+}(0,\mu_j,\sigma_j) = 1 - \Phi(0,\mu_j,\sigma_j)
- **T3.** The columns of S^{\*} have pairwise correlations
  R^{\*}\_{j\_{1}j\_{2}} given by \Phi^{+} and the upper-tail bivariate
  Gaussian cumulative distribution function \Phi^{+}\_2 as follows:
  \begin{equation\*} R^{\*}\_{j\_{1}j\_{2}} = \frac{
  \Phi^{+}\_{2}(\langle
  0,0\rangle,\langle\mu\_{j_1},\mu\_{j_2}\rangle,\Sigma\_{\langle
  j\_{1}, j\_{2}\rangle}) - \mu^{\*}\_{j_1}\mu^{\*}\_{j_2} }{
  \sqrt{\mu^{\*}\_{j_1}-{\mu^{\*}\_{j_1}}^2}
  \sqrt{\mu^{\*}\_{j_2}-{\mu^{\*}\_{j_2}}^2} } \end{equation\*}

## Proof (T1)

Let \Sigma^\prime=DR and note that by the definition of matrix
multiplication: \begin{equation\*} \Sigma^\prime\_{ij} =
\sum\_{k=1}^{m}D\_{ik}R\_{kj} \end{equation\*} By (D2), D\_{ik} = 0
except when i=k, in which case D\_{ik}=\sigma_i. Hence,
\Sigma^{\prime}\_{ij}=\sigma\_{i}R\_{ij}. Note that
\Sigma=\Sigma^{\prime}D, and thus: \begin{equation\*} \Sigma\_{ij} =
\sum\_{k=1}^{m}\Sigma^{\prime}\_{ik}D\_{kj} \end{equation\*} As before,
D\_{kj}=0 except when k=j, in which case D\_{kj}=\sigma\_{j}. Hence,
\Sigma\_{ij}=\Sigma^{\prime}\_{ij}\sigma\_{j}=\sigma\_{i}R\_{ij}\sigma\_{j}.
From the identity \mathrm{corr}(X,Y)\sigma_X\sigma_Y=\mathrm{cov}(X,Y)
and (S3), follows that \Sigma\_{ij}=\mathrm{cov}(s_i,s_j). Thus, \Sigma
is the covariance matrix for the columns of S. \blacksquare

## Proof (T2)

By the definition of \Phi^{+}, for 1\leq i\leq n: \begin{equation\*}
\Phi^{+}(0,\mu_j,\sigma_j)=\mathrm{P}(s\_{ij}\>
0\\\|\\s\_{ij}\sim\mathcal{N}(\mu_j,\sigma_j)) \end{equation\*} Hence by
(D3), \Phi^{+}(0,\mu_j,\sigma_j) is the number of values s^{\*}\_{ij}
which equal one, over n. This value is, by definition (S2, D1), the mean
\mu_j^{\*}. \blacksquare

## Proof (T3)

Consider again the identity: \begin{equation\*}
\mathrm{corr}(X,Y)=\frac{\mathrm{cov}(X,Y)}{\sigma\_{X}\sigma\_{Y}}
\end{equation\*} This identity implies that: \begin{equation\*}
R^{\*}\_{j\_{1}j\_{2}} = \frac{
\mathrm{cov}(s^{\*}\_{j\_{1}},s^{\*}\_{j\_{2}}) }{
\sigma^{\*}\_{j\_{1}}\sigma^{\*}\_{j\_{2}} } \end{equation\*} Given that
s^{\*}\_j\in\\0,1\\^{m} (D3), it follows that
s^{\*}\_{ij}={s^{\*}\_{ij}}^2 for all 1\leq i\leq n, so
\mu({s^{\*}\_j}^2)=\mu(s^{\*}\_j). Hence, given the identity:
\sigma(X)=\sqrt{\mu(X^2)-{\mu(X)}^2} it follows that: \begin{equation\*}
\sigma^{\*}\_{j} = \sqrt{\mu^{\*}\_{j}-{\mu^{\*}\_{j}}^2}
\end{equation\*} Therefore: \begin{equation\*} R^{\*}\_{j\_{1}j\_{2}} =
\frac{ \mathrm{cov}(s^{\*}\_{j\_{1}},s^{\*}\_{j\_{2}}) }{
\sqrt{\mu^{\*}\_{j_1}-{\mu^{\*}\_{j_1}}^2}
\sqrt{\mu^{\*}\_{j_2}-{\mu^{\*}\_{j_2}}^2} } \end{equation\*} Given (T2)
and the identity \mathrm{cov}(X,Y)=\mu(XY)-\mu(X)\mu(Y), it follows that
\mathrm{cov}(s^{\*}\_{j\_{1}},s^{\*}\_{j\_{2}}) =
\mu(s^{\*}\_{j\_{1}}s^{\*}\_{j\_{2}}) - \mu^{\*}\_{j_1}\mu^{\*}\_{j_2}.
Similar to proof of (T2), by the definition of \Phi^{+}, for 1\leq i\leq
n:\
\begin{gather\*} \Phi^{+}\_{2}(\langle
0,0\rangle,\langle\mu\_{j_1},\mu\_{j_2}\rangle,\Sigma\_{\langle j\_{1},
j\_{2}\rangle}) = \\ \mathrm{P}(s\_{ij\_{1}}\> 0, s\_{ij\_{2}} \> 0
\\\|\\ s\_{ij\_{1}}, s\_{ij\_{2}} \sim
\mathrm{MVN}(\langle\mu\_{j_1},\mu\_{j_2}\rangle,\Sigma\_{\langle
j\_{1}, j\_{2}\rangle})) \end{gather\*} Hence by (D3),
\Phi^{+}\_{2}(\langle
0,0\rangle,\langle\mu\_{j_1},\mu\_{j_2}\rangle,\Sigma\_{\langle j\_{1},
j\_{2}\rangle}) gives the number of rows i such that
s\_{ij\_{1}}=s\_{ij\_{2}}=1, over n. This value is equal to
\mu(s^{\*}\_{j\_{1}}s^{\*}\_{j\_{2}}). \blacksquare
