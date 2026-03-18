from mpi4py import MPI
import numpy as np

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

print("rank ", rank)

if rank == 0:
    sum = comm.recv(source=2)
    print(sum)
else:
    sum = comm.send(rank, dest=0)


