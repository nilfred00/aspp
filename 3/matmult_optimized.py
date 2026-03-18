# Program to multiply two matrices using numpy
import random
import numpy as np

@profile
def main():

    N = 250

    # NxN matrix
    X = np.random.rand(N, N)*100

    # Nx(N+1) matrix
    Y = np.random.rand(N, N+1)*100

    # multiply the two matrices 
    result = np.matmul(X, Y)

    # print the results
    print(result)

main()