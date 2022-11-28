import os, json
from collections import defaultdict

RESULTS_DIR = "results"
OUTPUT_FILE = os.path.join(RESULTS_DIR, "summary.json")

runs = []
for run in os.listdir(RESULTS_DIR):
    run_path = os.path.join(RESULTS_DIR, run)
    if os.path.isdir(run_path):
        runs.append((run, run_path))

print(f"Summarizing {len(runs)} runs...")

data = []
for run, run_path in runs:
    print(f"  Summarizing '{run}'")
    run_data = {}

    # Load the run.json file and ut it in the run_data dict.
    run_json_path = os.path.join(run_path, "run.json")
    if os.path.isfile(run_json_path):
        with open(run_json_path) as json_file:
            d = json.load(json_file)
            run_data.update(d)
    else:
        print("  WARNING: no run.json file in this directory.")
        continue

    groups = {}
    instance_group = {}

    # Scan the directory for "group_<groupname>.json" files.
    for file in os.listdir(run_path):
        if file.startswith("group_") and file.endswith(".json"):
            group_name = file[len("group_") :][: -len(".json")]
            print(f"  Found group name '{group_name}'.")
            group_json_path = os.path.join(run_path, file)
            with open(group_json_path) as json_file:
                d = json.load(json_file)
                groups[group_name] = {}
                for instance_obj in d["instances"]:
                    t_terminate = instance_obj["t_terminate"] if "t_terminate" in instance_obj else "?"
                    lb = instance_obj["solver_result"]["lb"] if "solver_result" in instance_obj and instance_obj["solver_result"] is not None else "?"
                    ub = instance_obj["solver_result"]["ub"] if "solver_result" in instance_obj and instance_obj["solver_result"] is not None else "?"
                    instance_name = instance_obj["name"]
                    if instance_name in instance_group:
                        print(
                            f"  WARNING: instance in multiple groups: '{instance_name}'"
                        )
                    instance_group[instance_name] = group_name
                    groups[group_name][instance_name] = {"solutions": [], "t_terminate": t_terminate, "lb": lb, "ub": ub}

    # Scan the directory for ".sol" files and parse their names for instance results.
    for file in os.listdir(run_path):
        if file.endswith(".sol"):
            filename_data = file[:-4]
            t, obj, instance_name = None, None, None
            split_filename = filename_data.split("_", maxsplit=2)

            t = float(split_filename[0])
            obj = float(split_filename[1])
            instance_name = split_filename[2]

            if not instance_name in instance_group:
                print(f"  WARNING: instance is not in any group: '{instance_name}'")
                continue

            instance_data = groups[instance_group[instance_name]][instance_name]
            instance_data["solutions"].append({"time": t, "obj": obj})

    n_solutions, n_instances = 0, 0
    for group in groups.values():
        for instance_name in group.values():
            n_instances += 1
            n_solutions += len(instance_name["solutions"])
            instance_name["solutions"].sort(key=lambda s: s["time"])

    print(
        f"  Found {n_solutions} solutions to {n_instances} instances in {len(groups)} groups."
    )
    run_data["groups"] = groups
    data.append(run_data)

with open(OUTPUT_FILE, "w", encoding="utf-8") as f:
    json.dump(data, f, ensure_ascii=False, indent=4)

# Also make a .js version so that we can use it from a local .html file.
with open(OUTPUT_FILE+ ".js", "w", encoding="utf-8") as f:
    f.write("var mipcomp_data = ")
    json.dump(data, f, ensure_ascii=False, indent=4)
    f.write(";\n")

print(f"Wrote summary file '{OUTPUT_FILE}'. See 'report.html'.")
