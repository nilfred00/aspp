from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

print(rank)
total = comm.reduce(rank, op=MPI.SUM, root=0)

if rank == 0:
    print(f"Total sum across all ranks: {total}")