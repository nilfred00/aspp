from mpi4py import MPI
import random

phrases = [
    "stay hydrated",
    "write cleaner code",
    "take a break",
    "optimize that loop",
    "commit your changes",
    "read the docs",
    "avoid premature optimization",
    "test your functions",
    "push to production carefully",
    "refactor that mess",
    "add more logging",
    "fix that bug",
    "keep it simple",
    "trust the process",
    "document everything",
    "never skip edge cases"
]

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

print(f"This is rank: {rank} telling you to {random.choice(phrases)}.")

