# svd_deconvolution
deconvolution of optical signal with singular value decomposition

An optical signal obtained from the discrete convolution between a tranfer function $\vec{h}$ wrapped in a Toeplitz matrix $A$ and an input signal $\vec{s}$ is given by

$$
\vec{m} =  A \vec{s}.
$$

The input signal is modified by the medium/environment, and this reaction is described by the transfer function $\vec{h}$.

Since $A$ isn't square it isnt inversible, we use a singular value decomposition instead to compute the pseudo inverse.
