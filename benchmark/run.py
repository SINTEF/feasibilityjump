import os, signal
import sys
import json
import time
import argparse
import subprocess
import platform

INSTANCES_DIR = "instances"
RESULTS_DIR = "results"
SOLVERS_DIR = "solvers"

# Arguments parser

parser = argparse.ArgumentParser(
    description="SINTEF MIP2022 competition perf test runner."
)
parser.add_argument("solver", help="The name or filename of the solver.")
parser.add_argument(
    "runname",
    help="The name of this test run. A directory with this name will be created to contain all output files.",
)
parser.add_argument(
    "--timeout",
    metavar="T",
    type=int,
    default=600,
    help="The timeout for each instance, in seconds. Defaults to 600.",
)
parser.add_argument("--all", action="store_true",
                    help="Run all instance groups.")
parser.add_argument("group", nargs="*", help="Group of instances to test.")
args = parser.parse_args()

# Make list of group names, or search instances directory for all group names if the '--all' argument is given.
groups = args.group
if args.all:
    groups = []
    for dir in os.listdir(INSTANCES_DIR):
        path = os.path.join(INSTANCES_DIR, dir)
        if os.path.isdir(path):
            groups.append(dir)

if len(groups) == 0:
    raise Exception("No instance groups!")

# Check that the output directory does not already exist, and then create it.
runname = args.runname

output_dir = os.path.join(RESULTS_DIR, runname)
try:
    os.makedirs(output_dir, exist_ok=False)
except FileExistsError:
    raise Exception("Run directory already exists. Delete it.")


# Create the solver function.
solver_exe = None
if os.path.isfile(args.solver):
    solver_exe = args.solver
elif os.path.isdir(os.path.join(SOLVERS_DIR, args.solver)):
    # Is this the directory name of a solver?
    # Try to find a python file or exe file inside
    for instance_name in [
        "main.py",
        "solver.py",
        "solver",
        args.solver,
        args.solver + ".exe",
        args.solver + ".py",
    ]:
        path = os.path.join(SOLVERS_DIR, args.solver, instance_name)
        if os.path.isfile(path):
            solver_exe = path

if solver_exe is None:
    raise Exception(f"Could not find solver executable for '{args.solver}'")

# If the solver is a python file, we need to run the python executable.
solver_exe_args = [solver_exe]
_, exe_extension = os.path.splitext(solver_exe)
print(f"EXE EXTENSION {exe_extension}")
if exe_extension == ".py":
    solver_exe_args = [sys.executable, solver_exe]
    if "PYTHONPATH" in os.environ and len(os.environ["PYTHONPATH"]) > 0:
        os.environ["PYTHONPATH"] += os.pathsep
    else:
        os.environ["PYTHONPATH"] = ""
    os.environ["PYTHONPATH"] += os.path.join(
        os.path.dirname(os.path.realpath(__file__)), os.pardir, "lib"
    )
    print(f"Set pythonpath to {os.environ['PYTHONPATH']}")

print(f"Using solver executable {solver_exe_args}")
timeout_secs = args.timeout


def run_solver(instance_file: str):
    start = time.time()
    solver_result = None
    try:
        args = solver_exe_args + [instance_file, output_dir]
        proc = subprocess.Popen(args, preexec_fn=os.setsid, stdout=subprocess.PIPE)
        print(f"    Running {args} as PID {proc.pid}...")
        try:
            #proc.wait(timeout=timeout_secs)
            outs, errs = proc.communicate(timeout=timeout_secs)
            for line in outs.decode('utf-8').splitlines():
                if line.startswith("{\""):
                    solver_result = json.loads(line)
        except subprocess.TimeoutExpired as e:
            # This is normal, the solver might run indefinitely, and that's ok.
            print(f"Killing process {proc.pid}.")
            os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
    except subprocess.CalledProcessError as e:
        print(f"Solver failed on instance {instance_file}: {e}")
    end = time.time()
    return end - start, solver_result

# Save parameters for the run to JSON


def get_git_revision_hash() -> str:
    return subprocess.check_output(['git', 'rev-parse', 'HEAD']).decode('ascii').strip()


git_rev = None
try:
    git_rev = get_git_revision_hash()
except:
    pass


run_data = {
    "runname": args.runname,
    "time": time.time(),
    "output_dir": output_dir,
    "solver": args.solver,
    "solver_exe": solver_exe_args,
    "timeout": timeout_secs,
    "groups": groups,
    "uname": dict(platform.uname()._asdict()),
    "git_rev": git_rev,
}

json_file = os.path.join(output_dir, "run.json")
with open(json_file, "w", encoding="utf-8") as f:
    json.dump(run_data, f, ensure_ascii=False, indent=4)

for group in groups:
    print(f"Running instance group '{group}'...")
    instances = []
    for instance_name in os.listdir(os.path.join(INSTANCES_DIR, group)):
        filename = os.path.join(INSTANCES_DIR, group, instance_name)
        if os.path.isfile(filename):
            _, extension = os.path.splitext(filename)
            if extension.lower() == ".mps" or extension.lower == ".lp":
                instances.append((instance_name, filename))

    # Save group data to JSON, so the report script knows which instances were run.

    print(f"  Group contains {len(instances)} instances.")
    instance_data = []
    for instance_name, path in instances:
        print(f"  Solving '{instance_name}'...")
        dt, solver_result = run_solver(path)
        instance_data.append((dt, solver_result))

    group_data = {
        "instances": [{"name": name, "path": path, "t_terminate": time, "solver_result": solver_result } for (name, path), (time, solver_result) in zip(instances, instance_times)],
    }

    with open(os.path.join(output_dir, f"group_{group}.json"), "w", encoding="utf-8") as f:
        json.dump(group_data, f, ensure_ascii=False, indent=4)
