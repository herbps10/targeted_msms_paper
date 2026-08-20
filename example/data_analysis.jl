
\subsection{Targeted Minimum Loss-based Estimation}
The motivation behind targeted estimators can be seen by analyzing the properties of a plug-in estimator. Suppose we have an initial estimate $\hat{P}$ of the data-generating distribution $\P$. The von Mises expansion reveals that the error of a plug-in estimator $\psi(\P)$ based on $\P$ decomposes as
\begin{align}
    \psi(\hat{P}) - \psi(\P) = \E_{P}\left\{ \D_{\hat{P}}(Z) \right\} + \R(\P, \hat{P}). 
\end{align}
The first term on the right hand side represents the first-order bias of the plug-in estimator. The second term is the second-order bias of the plug-in estimator. To show consistency of the plug-in estimator, we would need to show that both terms converge to zero asymptotically, which is difficult to ensure in a general setting. 

One approach is to \textit{de-bias} the plug-in estimate by estimating the first-order bias via empirical mean of the estimated EIF and subtracting it from the initial plug-in estimate:
\begin{align}
    \hat{\psi}^{\mathsf{onestep}} \equiv \psi(\hat{P}) - \frac{1}{n}\sum_{i=1}^n \D_{\hat{P}}(Z_i).
\end{align}
If conditions can be found to control the second-order remainder term $\R(\P, \hat{P})$, it is possible to show consistency and asymptotic normality. Such estimators are typically referred to as \textit{one-step estimators} because they are analogous to a a one-step Newton-Raphson update, with the EIF playing the role of the gradient of the functional \citep{fisher2021}. Under appropriate rate conditions on the nuisance estimators, $\hat{\psi}^{\mathsf{onestep}}$ is asymptotically normal and efficient. A drawback of this approach is that estimator is not guaranteed to lie in the parameter space; in our case, for example, a one-step estimate of the random replacement parameter may fall outside of $[0, 1]$. 

TMLE takes a different approach by carefully \textit{fluctuating} the initial estimate $\hat{P}$ of $\P$ along a low-dimensional parametric submodel in such a way that the updated estimate $\hat{P}^*$ (approximately) solves the empirical EIF estimating equation; that is,
\begin{align}
    \frac{1}{n}\sum_{i=1}^n \D_{\hat{P}^*}(Z_i) \approx 0.
\end{align}
A targeted estimate of $\psi_{a}^{\mathsf{replacement}}$ is then formed by plugging $\hat{P}^*$ into the parameter $\psi_{a}^{\mathsf{replacement}}(\cdot)$. As a substitution estimator, the targeted estimate is guaranteed by construction to lie in the parameter space. Constructing an appropriate fluctuation scheme that zeros the empirical first-order term is primarily a mathematical exercise; we give a complete description of the proposed algorithm in Appendix~\ref{appendix:tmle}. In terms of consistency, the form of the EIF implies that the resulting estimator will be consistent if \textit{either} the outcome regression $\mu$ or the attempt assignment model $\pi$ is estimated consistently, a result commonly referred to as \textit{double robustness}. 

By design, the targeted estimate eliminates the leading (first-order) term in the von Mises expansion. The remaining task is to establish conditions under which the second-order remainder term $\R(\P, \hat{P}^*)$ is asymptotically negligible. In our setting, the form the remainder can be shown to have a product form, involving the estimation errors of the outcome regression $\mu$ and the attempt assignment model $\pi$. As a result, the targeted estimator can be shown to be asymptotically normal and efficient as long as $\hat{\mu}$ and $\hat{\pi}$ converge at rates faster than $n^{-1/4}$. This result, referred to as \textit{rate double-robustness}, is the crucial finding that allows the use of data-adaptive nuisance estimation methods that converge at slower than parametric rates.

Our estimation strategy relies crucially on \textit{cross-fitting}. Cross-fitting proceeds by  splitting the data into multiple disjoint folds. For each fold, nuisance estimators are fit using observations in the complementary folds (the training set), and out-of-sample predictions are produced for observations in the held-out (validation) fold. Iterating over folds yields a complete set of nuisance predictions, where each prediction is computed out of sample. Cross-fitting allows us to use arbitrary, potentially complex nuisance estimators. Without cross-fitting, asymptotic normality results typically require restricting the complexity of the nuisance estimators. A classical way to achieve this is by assuming nuisance estimators fall in a Donsker class, verified through entropy arguments \citep{vanderVaart1996}. Cross-fitting obviates these requirements, and is now  a standard practice in semi-parametric estimation.

We are now  ready to formally state the statistical properties of the cross-fitted, targeted estimator of $\psi_{a}^{\mathsf{replacement}}$.
\begin{theorem}[Consistency and asymptotic normality]
    Let $\| \cdot \|$ denote the $L_2(\P)$ norm, and $\hat{\psi}_a^{\mathsf{replacement}}$ be the TMLE estimate of $\psi_{a}^{\mathsf{replacement}}$. Suppose that either $\| \hat{\pi} - \pi \| = o_{P}(1)$ or $\| \hat{\mu} - \mu \| = o_{P}(1)$. Then 
    \begin{align}
        \hat{\psi}_{a}^{\mathsf{replacement}} - \psi_{a}^{\mathsf{replacement}} = o_{P}(1) .
    \end{align}
    Suppose that both $\| \hat{\pi} - \pi \| = o_{P}(n^{-1/4})$ and $\| \hat{\mu} - \mu \| = o_{P}(n^{-1/4})$. Then 
    \begin{align}
        \sqrt{n}\left( \hat{\psi}_a^{\mathsf{replacement}} - \psi_{a}^{\mathsf{replacement}} \right) \leadsto N\left(0, \E_{P}\left[ \D^{\mathsf{replacement}}_{P}(Z)^2 \right]\right).
    \end{align}
\end{theorem}