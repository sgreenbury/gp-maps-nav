"""Python utilities for scripts."""
import uuid


def get_unique_id(integer=False):
    """
    Use uuid1 or uuid4 from lib to generate unique ID
    https://docs.python.org/3/library/uuid.html
    """
    if not integer:
        return str(uuid.uuid4()).replace("-", "")
    # Return a seed that can be used as an unsigned long
    return str(uuid.uuid4().int)[:9]


def proc_com(procs):
    """
    Communicate with open processes in procs list, holding loop until
    process completes.
    """
    for proc in procs:
        proc.communicate()


def print_args(args):
    """Print args."""
    for arg, val in vars(args).items():
        try:
            print(f"{arg:>39} : {val:>39}")
        except Exception:
            pass
