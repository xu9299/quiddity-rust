from pathlib import Path
import sys
import run
import arm_resources_run as common

def orders(stages):
    assert tuple(stages) == common.STAGES
    return [list(row) for row in common.orders()]

run.orders = orders
run.measure(*(Path(x).resolve() for x in sys.argv[1:]))
