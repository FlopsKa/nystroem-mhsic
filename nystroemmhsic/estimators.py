import numpy as np
from scipy.linalg import sqrtm
from functools import reduce
from nystroemmhsic.kernels import RBF

import fsic.data as data
import fsic.kernel as kernel
import fsic.indtest as it

class Estimator:
    def hsic(self, X, Y, kernelX, kernelY):
        pass

class MWEstimator:
    def mw_hsic(Xs, kernels):
        pass


class HSIC(Estimator, MWEstimator):
    """Naive HSIC estimator.
    """
    def __init__(self, kernels):
        self.kernels = kernels
     

    def hsic(self, X, Y, ):
        kernelX = self.kernels[0]
        kernelY = self.kernels[1]
        n = len(X)
        K = kernelX.k(X)
        L = kernelY.k(Y)

        H = np.identity(n) - (1 / n)

        KH = np.matmul(K, H)
        LH = np.matmul(L, H)
        KHLH = np.matmul(KH, LH)

        return 1 / (n**2) * np.trace(KHLH)

    def mw_hsic(self, Xs):
        n = len(Xs[0])
        one_n = 1 / n * np.ones((n, 1))
        K_XXs = [kernel.k(X) for X,kernel in zip(Xs,self.kernels)]
        A = one_n.T @ reduce(np.multiply, K_XXs) @ one_n
        B = reduce(np.multiply, [one_n.T @ K_i @ one_n for K_i in K_XXs])
        C = one_n.T @ reduce(np.multiply, [K_i @ one_n for K_i in K_XXs])
        return (A + B - 2 * C)[0][0]


class LargeScaleHSICNy(Estimator):
    """Nyström approximation for HSIC with $M=2$ components.
    """
    def __init__(self, kernelX, kernelY, num_nystrom_samples, reg_lambda=1e-8, seed=1234):
        """
        Arguments:
        ----------
        - num_nystrom_samples_func: a function with one parameter `n` which returns the number of samples to use, e.g., `np.sqrt`.
        - reg_lambda: regularization parameter
        - seed: seed for the rng.
        """
        self.kernelX = kernelX
        self.kernelY = kernelY
        self.samples = num_nystrom_samples
        self.reg_lambda = reg_lambda
        self.rng = np.random.default_rng(seed=seed)

    def hsic(self, X, Y):
        n = len(X)
        n_prime = self.samples
        idx = self.rng.integers(n,size=n_prime)
        X_tilde = X[idx]
        Y_tilde = Y[idx]

        K_X_nn_prime = self.kernelX.k(X,X_tilde)
        K_Y_nn_prime = self.kernelY.k(Y,Y_tilde)
        K_X_n_primen_prime = self.kernelX.k(X_tilde)
        K_Y_n_primen_prime = self.kernelY.k(Y_tilde)

        one_n = np.ones((n,1))
        H = np.identity(n) - 1/n

        Phi_X = K_X_nn_prime @ np.real(self.inv(sqrtm(K_X_n_primen_prime)))
        Phi_Y = K_Y_nn_prime @ np.real(self.inv(sqrtm(K_Y_n_primen_prime)))
        
        return 1 / (n**2) * np.linalg.norm((H @ Phi_X).T @ (H @ Phi_Y))**2

    def inv(self, A):
        return np.linalg.inv(A + self.reg_lambda * np.identity(len(A)))

    
class LargeScaleHSICRFF(Estimator):
    def __init__(self, kernelX, kernelY, num_rff_samples, seed=1234):
        self.kernelX = kernelX
        self.kernelY = kernelY
        self.samples = num_rff_samples
        self.rng = np.random.default_rng(seed=seed)

    def hsic(self, X, Y):
        n = len(X)
        s = self.samples
        H = np.identity(n) - 1/n
        
        X_rff = self.kernelX.Z(X, num_rff=s, rng=self.rng)
        Y_rff = self.kernelY.Z(Y, num_rff=s, rng=self.rng)

        return 1/(n**2) * np.linalg.norm(Y_rff.T @ H @ X_rff)**2
  
class HSICNy(Estimator, MWEstimator):
    def __init__(self, kernels, num_nystrom_samples, reg_lambda=1e-8, seed=1234):
        """
        Arguments:
        ----------
        - num_nystrom_samples_func: a function with one parameter `n` which returns the number of samples to use, e.g., `np.sqrt`.
        - reg_lambda: regularization parameter
        - seed: seed for the rng.
        """
        self.kernels = kernels
        self.reg_lambda = reg_lambda
        self.n_prime = num_nystrom_samples
        self.rng = np.random.default_rng(seed=seed)

    def hsic(self, X, Y):
        return self.mw_hsic(Xs = [X,Y]) 

    def mw_hsic(self, Xs):
        n = len(Xs[0])
        n_prime = self.n_prime
        if n_prime > n:
            return 0
        #idx = self.rng.integers(n,size=n_prime)
        idx = self.rng.choice(np.arange(n), size=n_prime, replace=True)
     
        # compute alpha for product kernel
        one_n = np.ones((n,1))
        mms = [kernel.k(X[idx]) for X, kernel in zip(Xs,self.kernels)] ## small Gram matrix
        mns = [kernel.k(X[idx], X) for X, kernel in zip(Xs,self.kernels)] ## rectangular matrix with Nyström samples
        mm = reduce(np.multiply,mms)
        mn = reduce(np.multiply,mns)
        alpha = 1/n * self.inv(mm) @ mn @ np.ones((n,1))
        # compute alpha for individual kernels
        alphas = [1/n * self.inv(mm) @ mn @ np.ones((n,1)) for mm, mn in zip(mms, mns)]
        
        A = alpha.T @ mm @ alpha
        B = reduce(np.multiply,[a.T @ mm @ a for a, mm in zip(alphas,mms)])
        C = alpha.T @ reduce(np.multiply,[mm @ a for a, mm in zip(alphas,mms)])
        
        return (A + B - 2 * C)[0][0]

    def inv(self, A):
        #return np.linalg.pinv(A, rcond=1e-15, hermitian=True)
        return np.linalg.inv(A + self.reg_lambda * np.identity(len(A)))

class FSICAdapter:
    def __init__(self, J, kernelX, kernelY, seed=1234,alpha=0.05):
        """
        Arguments:
        ----------
        - J: number of test point locations
        - seed: seed for the rng.
        """
        self.kernelX = kernelX
        self.kernelY = kernelY
        self.J = J
        self.seed = seed
        self.alpha = alpha
        gamma_x, gamma_y = self.kernelX.gamma, self.kernelY.gamma
        self.k = kernel.KGauss(gamma_x)
        self.l = kernel.KGauss(gamma_y)
        

    def hsic(self, X, Y):
        pdata = data.PairedData(X,Y)
        tr, te = pdata.split_tr_te(tr_proportion=0.5, seed=self.seed)
        J = self.J
        

        op = {
            'n_test_locs':J,  # number of test locations
            'max_iter':200, # maximum number of gradient ascent iterations
            'V_step':1, # Step size for the test locations of X
            'W_step':1, # Step size for the test locations of Y
            'gwidthx_step':1, # Step size for the Gaussian width of X
            'gwidthy_step':1, # Step size for the Gaussian width of Y
            'tol_fun':1e-4, # Stop if the objective function does not increase more than this
            'seed':self.seed+7 # random seed
        }
        op_V, op_W, op_gwx, op_gwy, info = it.GaussNFSIC.optimize_locs_widths(tr, self.alpha, **op )
        nfsic_opt = it.GaussNFSIC(op_gwx, op_gwy, op_V, op_W, self.alpha)
        
        
    
        return nfsic_opt.compute_stat(te)

    def test(self, X,Y, permutations=250,alpha=.05):
        pdata = data.PairedData(X,Y)
        tr, te = pdata.split_tr_te(tr_proportion=0.5, seed=self.seed)
        J = self.J
        op = {
            'n_test_locs':J,  # number of test locations
            'max_iter':100, # maximum number of gradient ascent iterations
            'V_step':1, # Step size for the test locations of X
            'W_step':1, # Step size for the test locations of Y
            'gwidthx_step':1, # Step size for the Gaussian width of X
            'gwidthy_step':1, # Step size for the Gaussian width of Y
            'tol_fun':1e-4, # Stop if the objective function does not increase more than this
            'seed':self.seed+7 # random seed
        }
        op_V, op_W, op_gwx, op_gwy, info = it.GaussNFSIC.optimize_locs_widths(tr, alpha, **op )
        nfsic_opt = it.GaussNFSIC(op_gwx, op_gwy, op_V, op_W, alpha, n_permute=permutations)
        return nfsic_opt.perform_test(te)["h0_rejected"]

    
class FSICAdapterNonOpt:
    def __init__(self, J, kernelX, kernelY, seed=1234,):
        """
        Arguments:
        ----------
        - J: number of test point locations
        - seed: seed for the rng.
        """
        self.kernelX = kernelX
        self.kernelY = kernelY
        self.J = J
        self.seed = seed
        gamma_x, gamma_y = self.kernelX.gamma, self.kernelY.gamma
        self.k = kernel.KGauss(gamma_x)
        self.l = kernel.KGauss(gamma_y)
        

    def hsic(self, X, Y):
        pdata = data.PairedData(X,Y)
        J = self.J
        V, W = it.GaussNFSIC.init_locs_marginals_subset(pdata, J, seed=self.seed)

        nfsic_med = it.NFSIC(self.k, self.l, V, W, reg='auto', n_permute=0, seed=self.seed)
        return nfsic_med.perform_test(pdata)["test_stat"]
